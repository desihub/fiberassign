#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "modules/htmTree.h"
#include "modules/kdTree.h"
#include "misc.h"
#include "feat.h"
#include "structs.h"
#include "global.h"

// t for global, time for local
Time t;

// Collecting information from input
// ------------------------------------------------------------
void collect_galaxies_for_all (const MTL & M,
                               const htmTree <struct target> & T, Plates & P,
                               const FP & pp,
                               const Feat & F) {
    // provides list of galaxies available to fiber k on tile j:
    // P[j].av_gals[k]
    init_time(t, "# Begin collecting available galaxies");
    double rad = F.PlateRadius * M_PI / 180.;

    #pragma omp parallel
    {
        #pragma omp master
        {
            printf(" ");
        }
        // Collects for each plate
        // start at jj=0 not id
        #pragma omp for
        for (int j = 0; j < F.Nplate; j++) {
            plate p = P[j];
            // Takes neighboring galaxies that fall on this plate
            std::vector <int> nbr = T.near(M, p.nhat, rad);
            // Projects thoses galaxies on the focal plane
            Onplates O;
            for (size_t gg = 0; gg < nbr.size(); gg++) {
                int g = nbr[gg];
                // struct onplate op = change_coords(M[g],p);
                struct onplate op = radec2xy(M[g], p);
                op.id = g;
                // Check that the target corresponds to the right program
                if (M[g].obsconditions & p.obsconditions) {
                    O.push_back(op);
                }
            }
            // Build 2D KD tree of those galaxies
            KDtree <struct onplate> kdT(O, 2);
            // For each fiber, finds all reachable galaxies within patrol
            // radius, thanks to the tree
            for (int k = 0; k < F.Nfiber; k++) {
                dpair X = pp[k].coords;
                std::vector <double> posit;
                posit.push_back(pp[k].fp_x);
                posit.push_back(pp[k].fp_y);
                std::vector <int> gals = kdT.near(&posit[0], 0.0, F.PatrolRad);
                for (size_t g = 0; g < gals.size(); g++) {
                    dpair Xg = projection(gals[g], j, M, P);
                    if (sq(Xg, X) < sq(F.PatrolRad) ) {
                        P[j].av_gals[k].push_back(gals[g]);
                        int q = pp[k].spectrom;
                        // better to make SS & SF assigned to fibers
                        if (M[gals[g]].SS) {
                            P[j].SS_av_gal[q].push_back(gals[g]);
                            P[j].SS_av_gal_fiber[k].push_back(gals[g]);
                        }
                        if (M[gals[g]].SF) {
                            P[j].SF_av_gal[q].push_back(gals[g]);
                            P[j].SF_av_gal_fiber[k].push_back(gals[g]);
                        }
                    }
                }
            }
        }
    }
    print_time(t, "# ... took :");
}

bool pairCompare (const std::pair <double, int> & firstElem,
                  const std::pair <double, int> & secondElem) {
    // used to sort galaxies by subpriority
    return firstElem.first < secondElem.first;
}

std::vector <int> sort_by_subpriority (MTL & M, std::vector <int> init) {
    // sorts list of galaxies by subpriority
    std::vector <int> out;
    std::vector <std::pair <double, int> > pairs;
    for (size_t gg = 0; gg < init.size(); ++gg) {
        int g = init[gg];
        std::pair <double, int> this_pair (M[g].subpriority, g);
        pairs.push_back(this_pair);
    }
    std::sort(pairs.begin(), pairs.end(), pairCompare);
    for (size_t gg = 0; gg < init.size(); ++gg) {
        out.push_back(pairs[gg].second);
    }
    return out;
}

bool int_pairCompare (const std::pair <int, int> & firstElem,
                      const std::pair <int, int> & secondElem) {
    // used to sort fibers by fib_num
    return firstElem.first < secondElem.first;
}

void collect_available_tilefibers (MTL & M, const Plates & P, const Feat & F) {
    // M[i].av_tfs is list of tile-fiber pairs available to galaxy i
    Time t;
    init_time(t, "# Begin computing available tilefibers");
    for (int j = 0; j < F.Nplate; j++) {
        for (int k = 0; k < F.Nfiber; k++) {
            for (size_t m = 0; m < P[j].av_gals[k].size(); m++) {
                // i is the id of the mth galaxy available to tile j and fiber
                // k
                int i = P[j].av_gals[k][m];
                // list of tile-fibers available to galaxy i
                M[i].av_tfs.push_back(pair(j, k) );
            }
        }
    }
    print_time(t, "# ... took :");
    int count_outside = 0;
    for (size_t g = 0; g < M.size(); ++g) {
        if ( ( M[g].av_tfs.size() == 0) && M[g].SF) {
            ++count_outside;
        }
    }
    printf("galaxies outside footprint %d\n", count_outside);
}

// Assignment sub-functions
// ------------------------------------------------------------------
// Allow (j,k) to observe g ?
inline bool ok_assign_g_to_jk (int g, int j, int k, const Plates & P,
                               const MTL & M, const FP & pp, const Feat & F,
                               const Assignment & A) {
    // if fiber is stuck or broken do not assign anything
    if (pp[k].stuck || pp[k].broken) {
        return false;
    }
    if (F.Collision) {
        for (size_t i = 0; i < pp[k].N.size(); i++) {
            if (g == A.TF[j][pp[k].N[i]]) {
                return false;
            }
        }
    }
    // Avoid 2 neighboring fibers observe the same galaxy (can happen only when
    // Collision=true)
    if (A.find_collision(j, k, g, pp, M, P, F) != -1) {
        return false;
    }
    // No collision
    return true;
    // doesn't require that jk is unassigned//doesn't require that g isn't
    // assigned already on this plate
    // use is_assigned_jg for this
}

// makes sure we don't exceed limit on SS and SF
inline bool ok_for_limit_SS_SF (int g, int j, int k, const MTL & M,
                                const Plates & P, const FP & pp,
                                const Feat & F) {
    bool is_SF = M[g].SF;
    bool too_many_SF = P[j].SF_in_petal[pp[k].spectrom] > F.MaxSF - 1;
    bool is_SS = M[g].SS;
    bool too_many_SS = P[j].SS_in_petal[pp[k].spectrom] > F.MaxSS - 1;
    return !(is_SF && too_many_SF) && !(is_SS && too_many_SS);
}

// Find, for (j,k), find the best galaxy it can reach among the possible ones
// Null list means you can take all possible kinds, otherwise you can only
// take, for the galaxy, a kind among this list
// Not allowed to take the galaxy of id no_g
inline int find_best (int j, int k, const MTL & M, const Plates & P,
                      const FP & pp, const Feat & F,
                      const Assignment & A) {
    int best = -1;
    int mbest = -1;
    int pbest = 0;
    double subpbest = 0.;
    List av_gals = P[j].av_gals[k];
    // For all available galaxies
    for (size_t gg = 0; gg < av_gals.size(); gg++) {
        int g = av_gals[gg];
        // if(ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){//don't assign SS, SF with
        // find_best 11/20/15
        if (!M[g].SS && !M[g].SF) {
            // Check whether it needs further observation
            int m = M[g].nobs_remain;
            if (m >= 1) {
                int prio = M[g].t_priority;
                double subprio = M[g].subpriority;
                // Takes it if better priority, or if same, if it needs more
                // observations,
                // so shares observations if two QSOs are close
                // If still tied, use subpriority
                if ( (prio > pbest) || ( (prio == pbest) && (m > mbest) ) ||
                     ( ( prio == pbest) && ( m == mbest) &&
                     ( subprio > subpbest) ) ) {
                    // Check that g is not assigned yet on this plate, or on
                    // the InterPlate around, check with ok_to_assign
                    int isa = A.is_assigned_jg(j, g, M, F);
                    int ok = ok_assign_g_to_jk(g, j, k, P, M, pp, F, A);
                    if ( (isa == -1) && ok) {
                        best = g;
                        pbest = prio;
                        mbest = m;
                        subpbest = subprio;
                    }
                }
            }
        }
    }
    return best;
}

// Tries to assign the fiber (j,k)
inline int assign_fiber (int j, int k, MTL & M, Plates & P, const FP & pp,
                         const Feat & F,
                         Assignment & A) {
    if (A.is_assigned_tf(j, k) ) {
        return -1;
    }
    int best = find_best(j, k, M, P, pp, F, A);
    if (best != -1) {
        A.assign(j, k, best, M, P, pp);
    }
    return best;
}

inline int assign_galaxy (int g,  MTL & M, Plates & P, const FP & pp,
                          const Feat & F, Assignment & A,
                          int jstart) {
    // Tries to assign the galaxy g to one of the used plates after jstart
    // jstart runs possibly to F.Nplate
    int jb = -1;
    int kb = -1;
    int unusedb = -1;
    Plist av_tfs = M[g].av_tfs;
    // All the tile-fibers that can observe galaxy g
    for (size_t tfs = 0; tfs < av_tfs.size(); tfs++) {
        int j = av_tfs[tfs].f;
        int k = av_tfs[tfs].s;
        // Check if the assignment is possible, if ok, if the tf is not used
        // yet, and if the plate is in the list
        if ( (jstart < j) &&
             !A.is_assigned_tf(j, k) &&
             ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) &&
             ok_for_limit_SS_SF(g, j, k, M, P, pp, F) ) {
            // unused fibers on this petal
            int unused = A.unused[j][pp[k].spectrom];
            if (unusedb < unused) {
                jb = j;
                kb = k;
                // observe this galaxy by fiber on petal with most free fibefs
                unusedb = unused;
            }
        }
    }
    if (jb != -1) {
        A.assign(jb, kb, g, M, P, pp);
        return 1;
    } else {
        return 0;
    }
}

// Takes an unassigned fiber and tries to assign it with the "improve"
// technique described in the doc
// not used for SS or SF
inline int improve_fiber (int jused_begin, int jused, int k, MTL & M,
                          Plates & P, const FP & pp, const Feat & F,
                          Assignment & A,
                          int no_g = -1) {
    int j = A.suborder[jused];
    if (!A.is_assigned_tf(j, k) ) {
        // Unused tilefiber (j,k)
        // maybe doesn't allow SS or SF
        int g_try = assign_fiber(j, k, M, P, pp, F, A);
        if ( (g_try != -1) && (g_try != -2) ) {
            return g_try;
        } else {
            // Improve
            int gb = -1;
            int bb = -1;
            int jpb = -1;
            int kpb = -1;
            int mb = -1;
            int pb = 1e3;
            int unusedb = -1;
            std::vector <int> av_g_init = P[j].av_gals[k];
            std::vector <int> av_g;
            // For all available galaxies within reach that are already
            // observed
            av_g = sort_by_subpriority(M, av_g_init);
            for (size_t i = 0; i < av_g.size(); i++) {
                // a galaxy accessible to j,k
                int g = av_g[i];
                if ( (g != -1) && (g != no_g) && !M[g].SS && !M[g].SF) {
                    // not SS or SF
                    if (ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                        // Which tile-fibers have taken g ?
                        // all tile-fibers that observe g in tiles from begin
                        // to end
                        Plist tfs =
                            A.chosen_tfs(g, F, A.suborder[jused_begin]);
                        for (size_t p = 0; p < tfs.size(); p++) {
                            // (jp,kp) currently assigned to galaxy g
                            int jp = tfs[p].f;
                            int kp = tfs[p].s;

                            // FIND BEST JP KP !!!
                            // best!=g because !A.assigned_pg(best)
                            int best = find_best(jp, kp, M, P, pp, F, A);
                            if ( (best != -1) &&
                                 ( (A.is_assigned_jg(j, g, M, F) == -1) ||
                                   ( jp == j) ) ) {
                                int prio = M[g].t_priority;
                                int m = M[g].nobs_remain;
                                // We take the most unused
                                int unused = A.unused[jp][pp[kp].spectrom];
                                if ( (prio > pb) ||
                                     ( (prio == pb) && (m > mb) ) ||
                                     ( ( prio == pb) && ( m == mb) &&
                                       ( unused > unusedb) ) ) {
                                    gb = g;
                                    bb = best;
                                    jpb = jp;
                                    kpb = kp;
                                    mb = m;
                                    pb = prio;
                                    unusedb = unused;
                                }
                            }
                        }
                    }
                }
            }
            // Modify assignment
            if (gb != -1) {
                A.unassign(jpb, kpb, gb, M, P, pp);
                A.assign(j, k, gb, M, P, pp);
                A.assign(jpb, kpb, bb, M, P, pp);
                return gb;
            }
        }
    }
    return -1;
}

// Assignment functions
// ---------------------------------------------------------
// Assign fibers naively
bool inverse_pairCompare (const std::pair <double, int> & firstElem,
                          const std::pair <double, int> & secondElem) {
    // used to sort galaxies by subpriority
    return firstElem.first > secondElem.first;
}

void simple_assign (MTL & M, Plates & P, const FP & pp, const Feat & F,
                    Assignment & A) {
    Time t;
    init_time(t, "# Begin simple assignment :");
    int countme = 0;
    for (int j = 0; j < F.Nplate; j++) {
        fprintf(stdout, "plate %d out of %d\n", j, F.Nplate);
        fflush(stdout);
        std::vector <std::pair <double, int> > pairs;
        for (int k = 0; k < F.Nfiber; k++) {

            List av_gals = P[j].av_gals[k];
            double Max_fiber_priority = 0.;
            // loop over the potential list to find pbest and subpbest
            // and then combine them to give the fiber_weight
            for (size_t gg = 0; gg < av_gals.size(); gg++) {
                //fprintf(stdout, "available gal %d out of %d\n", gg, av_gals.size());
                //    fflush(stdout);
                int g = av_gals[gg];
                // not SS and SF
                if (!M[g].SS && !M[g].SF) {
                    int m = M[g].nobs_remain;
                    // unobserved
                    if (m >= 1) {
                        int prio = M[g].t_priority;
                        double subprio = M[g].subpriority;
                        double fiber_priority = (double)prio + subprio;
                        if (fiber_priority > Max_fiber_priority) {
                            Max_fiber_priority = fiber_priority;
                        }
                    }
                }
            }
            std::pair <double, int> this_pair(Max_fiber_priority, k);
            pairs.push_back(this_pair);
        }
        std::sort(pairs.begin(), pairs.end(), inverse_pairCompare);
        std::vector <int> fiber_loop_order;
        for (int k = 0; k < F.Nfiber; k++) {
            fiber_loop_order.push_back(pairs[k].second);
        }
        int best = -1;
        for (int k = 0; k < F.Nfiber; k++) {
            // Fiber
            best = assign_fiber(j, fiber_loop_order[k], M, P, pp, F, A);
            if ( ( best != -1) || ( best != -2) ) {
                countme++;
            }
        }
    }
    print_time(t, "# ... took :");
    printf(" countme %d \n", countme);
}

void improve (MTL & M, Plates & P, const FP & pp, const Feat & F,
              Assignment & A, int jused_start) {
    // jstart is in list from 0 to F.NUsedplate
    Time t;
    init_time(t, "# Begin improve :");
    int improvements = 0;
    for (int jused = jused_start; jused < F.NUsedplate; jused++) {
        for (int k = 0; k < F.Nfiber; k++) {
            int worked = improve_fiber(jused_start, jused, k, M, P, pp, F, A);
            if (worked != -1) improvements++;
        }
    }
    printf(" improvements  %d\n", improvements);
    print_time(t, "# ... took :");
}

// If there are galaxies discovered as fake for example, they won't be observed
// several times in the plan
// has access to G,not just M, because it needs to know the truth
void new_replace (int j, int p, MTL & M, Plates & P, const FP & pp,
                  const Feat & F, Assignment & A) {
    // do standard stars,going through priority classes from least to most
    // keep track of reassignments
    int reassign_SS = 0;
    int reassign_SF = 0;
    int add_SS = 0;
    int add_SF = 0;
    // can get all available SS,SF on plate from P[j].av_gals_plate restricting
    // to plate p
    for (size_t c = 0;
        P[j].SS_in_petal[p] < F.MaxSS && c < M.priority_list.size(); ++c ) {
        // try to do this for lowest priority first
        // standard stars on this petal
        std::vector <int> gals_init = P[j].SS_av_gal[p];
        std::vector <int> gals;
        // sort these by priority
        gals = sort_by_subpriority(M, gals_init);
        for (size_t gg = 0; gg < gals.size(); ++gg) {
            // all the standard stars on this petal
            // a standard star
            int g = gals[gg];
            if (A.is_assigned_jg(j, g) == -1) {
                // not assigned on this tile
                // all tiles and fibers that reach g
                Plist tfs = M[g].av_tfs;
                // quit after we've used this SS
                int done = 0;
                std::vector <int> could_be_replaced;
                std::vector <int> sorted_could_be_replaced;
                for (size_t i = 0; i < tfs.size() && done == 0; ++i) {
                    // al the tile-fibers that can reach this standard star
                    if ( (tfs[i].f == j) && (pp[tfs[i].s].spectrom == p) ) {
                        // a tile fiber from this petal
                        // we know g can be reached by this petal of plate j
                        // and fiber k
                        int k = tfs[i].s;
                        // what is now at (j,k)  g_old can't be -1 or we would
                        // have used it already in assign_sf
                        int g_old = A.TF[j][k];
                        // require that this galaxy is used only once to keep
                        // things simple
                        if ( (g_old != -1) && !M[g_old].SS && !M[g_old].SF &&
                             ( A.GL[g_old].size() == 1) ) {
                            if ( ((size_t)M[g_old].priority_class == c) &&
                                ( A.is_assigned_jg(j, g, M, F) == -1) &&
                                ok_for_limit_SS_SF(g, j, k, M, P, pp, F) &&
                                ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                                // right priority; this SS not already assigned
                                // on this plate
                                could_be_replaced.push_back(g_old);
                            }
                        }
                    }
                }
                if (could_be_replaced.size() > 0) {
                    // now have list of possible targets to be replaced with
                    // lowest priority
                    // use subpriority to choose
                    std::vector <int> sorted_could_be_replaced =
                        sort_by_subpriority(M, could_be_replaced);
                    int g_chosen = sorted_could_be_replaced[0];
                    int j = A.GL[g_chosen][0].f;
                    int k = A.GL[g_chosen][0].s;
                    A.unassign(j, k, g_chosen, M, P, pp);
                    // try to assign, increment counter
                    reassign_SS += assign_galaxy(g_chosen, M, P, pp, F, A, j);
                    A.assign(j, k, g, M, P, pp);
                    add_SS++;
                    done = 1;
                }
            }
        }
    }
    for (size_t c = 0;
        P[j].SF_in_petal[p] < F.MaxSF && c < M.priority_list.size(); ++c ) {
        // try to do this for lowest priority
        // sky fibers on this petak
        std::vector <int> gals_init = P[j].SF_av_gal[p];
        // sort these by subpriority
        std::vector <std::pair <double, int> > galaxy_pairs;
        for (size_t gg = 0; gg < gals_init.size(); ++gg) {
            int g = gals_init[gg];
            std::pair <double, int> this_pair (M[g].subpriority, g);
            galaxy_pairs.push_back(this_pair);
        }
        std::sort(galaxy_pairs.begin(), galaxy_pairs.end(), pairCompare);
        std::vector <int> gals;
        for (size_t gg = 0; gg < gals_init.size(); ++gg) {
            gals.push_back(galaxy_pairs[gg].second);
        }
        for (size_t gg = 0; gg < gals.size(); ++gg) {
            int g = gals[gg];  // a sky fiber
            if (A.is_assigned_jg(j, g) == -1) {
                Plist tfs = M[g].av_tfs;
                int done = 0;
                std::vector <int> could_be_replaced;
                for (size_t i = 0; i < tfs.size() && done == 0; ++i) {
                    if ( (tfs[i].f == j) && (pp[tfs[i].s].spectrom == p) ) {
                        // g is accesible to j,k
                        // we know g can be reached by this petal of plate j
                        // and fiber k
                        int k = tfs[i].s;
                        // what is now at (j,k)
                        int g_old = A.TF[j][k];
                        if ( (g_old != -1) && !M[g_old].SS && !M[g_old].SF &&
                             ( A.GL[g_old].size() == 1) ) {
                            if ( ((size_t)M[g_old].priority_class == c) &&
                                ( A.is_assigned_jg(j, g, M, F) == -1) &&
                                ok_for_limit_SS_SF(g, j, k, M, P, pp, F) &&
                                ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                                could_be_replaced.push_back(g_old);
                            }
                        }
                    }
                }
                // now have list of possible targets to be replaced with lowest
                // priority
                // use subpriority to choose
                sort_by_subpriority(M, could_be_replaced);
                if (could_be_replaced.size() > 0) {
                    int g_chosen = could_be_replaced[0];
                    int j = A.GL[g_chosen][0].f;
                    int k = A.GL[g_chosen][0].s;
                    A.unassign(j, k, g_chosen, M, P, pp);
                    // try to assign
                    reassign_SF += assign_galaxy(g_chosen, M, P, pp, F, A, j);
                    A.assign(j, k, g, M, P, pp);
                    add_SF++;
                    done = 1;
                }
            }
        }
    }
}

void assign_unused (int j, MTL & M, Plates & P, const FP & pp, const Feat & F,
                    Assignment & A) {
    // Tries to assign remaining fibers in tile jth tile with galaxies on it
    // even taking objects observed later
    // js is a tile with galaxies on it
    for (int k = 0; k < F.Nfiber; k++) {
        if (!A.is_assigned_tf(j, k) ) {
            int best = -1;
            int mbest = -1;
            int pbest = 0;
            int jpb = -1;
            int kpb = -1;
            double subpbest = 0;
            // all available galaxies for this fiber k
            std::vector <int> av_gals = P[j].av_gals[k];
            std::vector <int> sorted_av_gals = sort_by_subpriority(M, av_gals);
            for (size_t gg = 0; gg < sorted_av_gals.size(); gg++) {
                // available galaxies
                int g = sorted_av_gals[gg];
                int m = M[g].nobs_remain;
                int prio = M[g].t_priority;
                double subprio = M[g].subpriority;
                if ( (prio > pbest) || ( (prio == pbest) && (m > mbest) ) ||
                     ( ( prio == pbest) && ( m == mbest) &&
                       ( subprio > subpbest) ) ) {
                    if ( (A.is_assigned_jg(j, g, M, F) == -1) &&
                         ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) &&
                         ok_for_limit_SS_SF(g, j, k, M, P, pp, F) ) {
                        // not assigned this plate or within excluded interval
                        for (size_t i = 0; i < A.GL[g].size(); i++) {
                            // GL[g].size() is number of tf that could observe
                            // g
                            int jp = A.GL[g][i].f;
                            int kp = A.GL[g][i].s;
                            if ( (j < jp) && (jpb < jp) ) {
                                // take best opportunity
                                best = g;
                                pbest = prio;
                                mbest = m;
                                subpbest = subprio;
                                jpb = jp;
                                kpb = kp;
                            }
                        }
                    }
                }
            }
            if (best != -1) {
                A.unassign(jpb, kpb, best, M, P, pp);
                A.assign(j, k, best, M, P, pp);
            }
        }
    }
}

void assign_sf_ss (int j, MTL & M, Plates & P, const FP & pp, const Feat & F,
                   Assignment & A) {
    // assigns sky fibers and standard stars to unused fibers
    for (int ppet = 0; ppet < F.Npetal; ppet++) {
        // for each petal
        int p = ppet;
        std::vector <int> SS_av_init = P[j].SS_av_gal[p];
        std::vector <int> SS_av;
        SS_av = sort_by_subpriority(M, SS_av_init);
        // use subpriority to order SS, SF
        std::vector <int> SF_av_init = P[j].SF_av_gal[p];
        std::vector <int> SF_av;
        SF_av = sort_by_subpriority(M, SF_av_init);
        if ( (SS_av.size() > 0) || (SF_av.size() > 0) ) {
            // look at fibers on this petal
            for (int kk = 0; kk < F.Nfbp; kk++) {
                int k = pp.fibers_of_sp[p][kk];
                // for a particular fiber on this petal consider SS and SF by
                // subpriority order
                std::vector <int> SS_av_k_init = P[j].SS_av_gal_fiber[k];
                std::vector <int> SF_av_k_init = P[j].SF_av_gal_fiber[k];
                std::vector <int> SS_av_k;
                std::vector <int> SF_av_k;
                SS_av_k = sort_by_subpriority(M, SS_av_k_init);
                SF_av_k = sort_by_subpriority(M, SF_av_k_init);
                if (A.TF[j][k] == -1) {
		   // this fiber isn't assigned yet
		   int done = 0;

		   for (size_t gg = 0; gg < SS_av_k.size() && done == 0; gg++) {
                        // SS on this with reach of this fiber
                       int g = SS_av_k[gg];
                        // if g isn't already assigned on this fiber and we
                        // still need it and it's ok to assign, then assign
                       if ( (A.is_assigned_jg(j, g, M, F) == -1) &&
                             ok_for_limit_SS_SF(g, j, k, M, P, pp, F) &&
                             ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                            A.assign(j, k, g, M, P, pp);
                            done = 1;
                       }
                    }
                    for (size_t gg = 0; gg < SF_av_k.size()
                        && done == 0; gg++) {
                        // SF on this petal
                        int g = SF_av_k[gg];
                        if ( (A.is_assigned_jg(j, g, M, F) == -1) &&
                             ok_for_limit_SS_SF(g, j, k, M, P, pp, F) &&
                             ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                            A.assign(j, k, g, M, P, pp);
                            done = 1;
                        }
                    }
                }
            }
            new_replace(j, p, M, P, pp, F, A);
        }
    }
}

void fill_unused_with_sf (int j, MTL & M, Plates & P, const FP & pp,
                          const Feat & F, Assignment & A) {
    // assigns sky fibers to fill all unused fibers
    for (int ppet = 0; ppet < F.Npetal; ppet++) {
        // for each petal
        int p = ppet;
        // use subpriority to order SF (Sky Fiber)
        std::vector <int> SF_av_init = P[j].SF_av_gal[p];
        std::vector <int> SF_av;
        SF_av = sort_by_subpriority(M, SF_av_init);
        if (SF_av.size() > 0) {
            // look at fibers on this petal
            for (int kk = 0; kk < F.Nfbp; kk++) {
                int k = pp.fibers_of_sp[p][kk];
                // for a particular fiber on this petal consider SF by
                // subpriority order
                std::vector <int> SF_av_k_init = P[j].SF_av_gal_fiber[k];
                std::vector <int> SF_av_k;
                SF_av_k = sort_by_subpriority(M, SF_av_k_init);
                if (A.TF[j][k] == -1) {
                    // this fiber isn't assigned yet
                    int done = 0;
                    for (size_t gg = 0; gg < SF_av_k.size()
                        && done == 0; gg++) {
                        // SF on this petal
                        int g = SF_av_k[gg];
                        if ( (A.is_assigned_jg(j, g, M, F) == -1) &&
                             ok_assign_g_to_jk(g, j, k, P, M, pp, F, A) ) {
                            A.assign(j, k, g, M, P, pp);
                            done = 1;
                        }
                    }
                }
            }
            new_replace(j, p, M, P, pp, F, A);
        }
    }
}

void redistribute_tf (MTL & M, Plates & P, const FP & pp, const Feat & F,
                      Assignment & A, int jused_start) {
    // diagnostic
    printf("start redistribute \n");
    Time t;
    init_time(t, "# Begin redistribute TF :");
    int red(0);
    // consider every occupied plate and every fiber
    Table Done = initTable(F.NUsedplate, F.Nfiber);
    for (int jused = jused_start; jused < F.NUsedplate; jused++) {
        int j = A.suborder[jused];
        for (int k = 0; k < F.Nfiber; k++) {
            if (Done[jused][k] == 0) {
                // current assignment of (j,k)  only look if assigned
                int g = A.TF[j][k];
                if ( (g != -1) && !M[g].SS && !M[g].SF) {
                    int jpb = -1;
                    int kpb = -1;
                    int unusedb = A.unused[j][pp[k].spectrom];
                    // all possible tile fibers for this galaxy
                    Plist av_tfs = M[g].av_tfs;
                    for (size_t i = 0; i < av_tfs.size(); i++) {
                        int jp = av_tfs[i].f;
                        int kp = av_tfs[i].s;
                        // unused for jp, spectrom[kp]
                        int unused = A.unused[jp][pp[kp].spectrom];
                        if (A.inv_order[jp] != -1) {
                            // necessary because underdense targets may leave
                            // some plates unused
                            if ( ( A.inv_order[jp] > F.NUsedplate) ||
                                 ( A.inv_order[jp] < 0) ) {
                                printf("**out range  %d\n", A.inv_order[jp]);
                            }
                            if (A.suborder[jused_start] <= jp) {
                                if (!A.is_assigned_tf(jp, kp) ) {
                                    if (Done[A.inv_order[jp]][kp] == 0) {
                                        if ( ok_assign_g_to_jk(g, jp, kp, P, M,
                                                               pp, F, A) ) {
                                            if (A.is_assigned_jg(jp, g, M, F)
                                                == -1) {
                                                if ( 0 < unused) {
                                                    if (unusedb < unused) {
                                                        // Takes the most
                                                        // unused petal
                                                        jpb = jp;
                                                        kpb = kp;
                                                        unusedb = unused;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (jpb != -1) {
                        A.unassign(j, k, g, M, P, pp);
                        A.assign(jpb, kpb, g, M, P, pp);
                        Done[A.inv_order[j]][k] = 1;
                        Done[A.inv_order[jpb]][kpb] = 1;
                        red++;
                    }
                }
            }
        }
    }
    printf("  %s redistributions of tile-fibers \n", f(red).c_str() );
    print_time(t, "# ... took :");
}

void fa_write (int j, str outdir, const MTL & M, const Plates & P,
               const FP & pp, const Feat & F, const Assignment & A) {
    // generate a quiet NaN to use for invalid entries.  We cannot
    // guarantee that we have C++11, so we can't use the nice functions
    // included in that standard...
    // const unsigned long maxU = ~0;
    // const double qNan = *((double*)&maxU);
    // double qNan = nan("");
    double testra, testdec;
    double arcsec = 1.0 / 3600.0;
    // constants for the filename length and fixed object
    // type length
    size_t cfilesize = 512;
    //size_t objtypelen = 8;
    size_t bricklen = 8;
    // check if the file exists, and if so, throw an exception
    char filename[cfilesize];
    int ret = snprintf(filename, cfilesize, "%s/tile_%05d.fits",
                       outdir.c_str(), P[j].tileid);
    struct stat filestat;
    ret = ::stat(filename, &filestat );
    if (ret == 0) {
        std::cerr << "ERROR: output file " << filename << " already exists"
            << std::endl;
        std::exit(1);
        // std::ostringstream o;
        // o << "output file " << filename << " already exists";
        // throw std::runtime_error(o.str().c_str());
    }
    // create the file
    int status = 0;
    fitsfile * fptr;
    fits_create_file (&fptr, filename, &status);
    fits_report_error (stderr, status);
    // Set up the schema for the table.  We explicitly malloc these
    // string arrays, since the CFITSIO API requires non-const pointers
    // to them (i.e. arrays of literals won't work).
    size_t ncols = 14;
    char * * ttype;
    char * * tform;
    char * * tunit;
    ttype = (char * *) malloc ( ncols * sizeof(char *) );
    tform = (char * *) malloc ( ncols * sizeof(char *) );
    tunit = (char * *) malloc ( ncols * sizeof(char *) );
    if ( !( ttype && tform && tunit ) ) {
        std::ostringstream o;
        o << "cannot allocate column info for binary table";
        throw std::runtime_error(o.str().c_str() );
    }
    for ( size_t c = 0; c < ncols; ++c ) {
        ttype[c] = (char *) malloc ( FLEN_VALUE * sizeof(char) );
        tform[c] = (char *) malloc ( FLEN_VALUE * sizeof(char) );
        tunit[c] = (char *) malloc ( FLEN_VALUE * sizeof(char) );
        if ( !( ttype[c] && tform[c] && tunit[c] ) ) {
            std::ostringstream o;
            o << "cannot allocate column info for binary table";
            throw std::runtime_error(o.str().c_str() );
        }
    }
    strcpy(ttype[0], "FIBER");
    strcpy(tform[0], "J");
    strcpy(tunit[0], "");
    strcpy(ttype[1], "LOCATION");
    strcpy(tform[1], "J");
    strcpy(tunit[1], "");
    strcpy(ttype[2], "NUMTARGET");
    strcpy(tform[2], "I");  // int not long
    strcpy(tunit[2], "");
    strcpy(ttype[3], "PRIORITY");
    strcpy(tform[3], "J");
    strcpy(tunit[3], "");
    strcpy(ttype[4], "TARGETID");
    strcpy(tform[4], "K");
    strcpy(tunit[4], "");
    strcpy(ttype[5], "DESI_TARGET");
    strcpy(tform[5], "K");
    strcpy(tunit[5], "");
    strcpy(ttype[6], "BGS_TARGET");
    strcpy(tform[6], "K");
    strcpy(tunit[6], "");
    strcpy(ttype[7], "MWS_TARGET");
    strcpy(tform[7], "K");
    strcpy(tunit[7], "");
    strcpy(ttype[8], "RA");
    strcpy(tform[8], "D");
    strcpy(tunit[8], "deg");
    strcpy(ttype[9], "DEC");
    strcpy(tform[9], "D");
    strcpy(tunit[9], "deg");
    strcpy(ttype[10], "XFOCAL_DESIGN");
    strcpy(tform[10], "E");
    strcpy(tunit[10], "mm");
    strcpy(ttype[11], "YFOCAL_DESIGN");
    strcpy(tform[11], "E");
    strcpy(tunit[11], "mm");
    strcpy(ttype[12], "BRICKNAME");
    snprintf(tform[12], FLEN_VALUE, "%dA", (int)bricklen);
    strcpy(tunit[12], "");
    strcpy(ttype[13], "FIBERMASK");
    strcpy(tform[13], "J");
    strcpy(tunit[13], "");
    char extname[FLEN_VALUE];
    strcpy(extname, "FIBER_ASSIGNMENTS");
    // create the table with the full size on disk.
    ret = fits_create_tbl(fptr, BINARY_TBL, F.Nfiber, ncols, ttype, tform,
                          tunit, extname, &status);
    fits_report_error(stderr, status);
    int tileid = P[j].tileid;
    double tilera = P[j].tilera;
    double tiledec = P[j].tiledec;
    fits_write_key(fptr, TINT, "TILEID", &(tileid), "Tile ID number", &status);
    fits_report_error(stderr, status);
    fits_write_key(fptr, TDOUBLE, "TILERA", &(tilera), "Tile RA [deg]",
                   &status);
    fits_report_error(stderr, status);
    fits_write_key(fptr, TDOUBLE, "TILEDEC", &(tiledec), "Tile DEC [deg]",
                   &status);
    fits_report_error(stderr, status);
    // get the number of rows to write for each internal FITS buffer.
    long optimal;
    ret = fits_get_rowsize(fptr, &optimal, &status);
    fits_report_error(stderr, status);
    // initialize arrays to the optimal number of rows for writing.
    //int tile_id[optimal];
    int fiber_id[optimal];
    int positioner_id[optimal];
    int num_target[optimal];
    int fibermask[optimal];
    //char objtype[optimal][objtypelen];
    char brickname[optimal][bricklen + 1];
    char * bn_tmp[optimal];
    //char * ot_tmp[optimal];
    for (int i = 0; i < optimal; i++) {
        //ot_tmp[i] = objtype[i];
        bn_tmp[i] = brickname[i];
    }
    long long target_id[optimal];
    long long desi_target[optimal];
    long long bgs_target[optimal];
    long long mws_target[optimal];
    double ra[optimal];
    double dec[optimal];
    double x_focal[optimal];
    double y_focal[optimal];
    // new
    int t_priority[optimal];
    std::vector <long long> potentialtargetid;
    long long offset = 0;
    long long n = optimal;
    while ( n == optimal ) {
        if ( offset + optimal > F.Nfiber ) {
            n = F.Nfiber - offset;
        }
        if ( n > 0 ) {
            for (long long i = 0; i < n; ++i) {
                long long fib = offset + i;
                int g = A.TF[j][fib];
                fiber_id[i] = fib;
                positioner_id[i] = pp[fib].location;
                num_target[i] = P[j].av_gals[fib].size();
                // target_id[i] = g; ********
                if (g > 0) {
                    target_id[i] = M[g].id;
                } else {
                    target_id[i] = -1;
                }
                if (g < 0) {
                    // For negative target ID, must still initialize
                    // memory buffer that will be written
                    t_priority[i] = 1;

                    if (pp[fib].stuck) {
                        fibermask[i] = FIBER_STUCK;
                        x_focal[i] = pp[fib].fp_x;
                        y_focal[i] = pp[fib].fp_y;
                        xy2radec(&(ra[i]), &(dec[i]), tilera, tiledec,
                                 x_focal[i], y_focal[i]);
                        desi_target[i] = 0;
                        bgs_target[i] = 0;
                        mws_target[i] = 0;
                        strncpy(brickname[i], "notbrick", bricklen + 1);
                    } else if (pp[fib].broken) {
                        fibermask[i] = FIBER_BROKEN;
                        x_focal[i] = pp[fib].fp_x;
                        y_focal[i] = pp[fib].fp_y;
                        xy2radec(&(ra[i]), &(dec[i]), tilera, tiledec,
                                 x_focal[i], y_focal[i]);
                        desi_target[i] = 0;
                        bgs_target[i] = 0;
                        mws_target[i] = 0;
                        strncpy(brickname[i], "notbrick", bricklen + 1);
                    } else {
                        fibermask[i] = FIBER_UNUSED;
                        x_focal[i] = pp[fib].fp_x;
                        y_focal[i] = pp[fib].fp_y;
                        xy2radec(&(ra[i]), &(dec[i]), tilera, tiledec,
                                 x_focal[i], y_focal[i]);
                        desi_target[i] = 0;
                        bgs_target[i] = 0;
                        mws_target[i] = 0;
                        strncpy(brickname[i], "notbrick", bricklen + 1);
                    }
                } else {
                    // we aren't supposed to know the kind  use priority
                    // instead
                    // strncpy(objtype[i], F.kind[G[g].id].c_str(),
                    // objtypelen);
                    fibermask[i] = FIBER_USED;
                    ra[i] = M[g].ra;
                    dec[i] = M[g].dec;
                    dpair proj = projection(g, j, M, P);
                    x_focal[i] = proj.f;
                    y_focal[i] = proj.s;
                    // testing that xy2radec is working with good precision.
                    xy2radec(&testra, &testdec, tilera, tiledec, x_focal[i],
                             y_focal[i]);
                    double dra =
                        (testra * cos(testdec * M_PI / 180.0) - ra[i] *
                            cos(dec[i] * M_PI / 180.0) ) / arcsec;
                    double ddec = (testdec - dec[i]) / arcsec;
                    if (fabs(dra) > 0.01) {
                        // arcsecond precision
                        fprintf(stderr,
                                "problem with xy2radec conversion [dRA (arcsec)]: %f\n",
                                dra);
                        fprintf(stderr, "[dDEC (arcsec)]: %f \n", ddec);
                        myexit(1);
                    }
                    if (fabs(ddec) > 0.01) {
                        // arcsecond precision
                        fprintf(stderr,
                                "problem with xy2radec conversion [dDEC]: %f\n",
                                ddec);
                        fprintf(stderr, "[dRA]: %f\n", dra);
                        myexit(1);
                    }
                    t_priority[i] = M[g].t_priority;  // new
                    desi_target[i] = M[g].desi_target;
                    bgs_target[i] = M[g].bgs_target;
                    mws_target[i] = M[g].mws_target;
                    strncpy(brickname[i], M[g].brickname, bricklen + 1);
                }
                // Store the potential targetids accesible to this fibre (the
                // actual targetid, not the index).
                for (size_t k = 0; k < P[j].av_gals[fib].size(); ++k) {
                    // MTL index for k'th target accessible to this fibre
                    int gal_idx = P[j].av_gals[fib][k];
                    if (gal_idx >= 0) {
                        potentialtargetid.push_back(M[gal_idx].id);
                    }
                }
            }
            fits_write_col(fptr, TINT, 1, offset + 1, 1, n, fiber_id, &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TINT, 2, offset + 1, 1, n, positioner_id,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TINT, 3, offset + 1, 1, n, num_target,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TINT, 4, offset + 1, 1, n, t_priority,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TLONGLONG, 5, offset + 1, 1, n, target_id,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TLONGLONG, 6, offset + 1, 1, n, desi_target,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TLONGLONG, 7, offset + 1, 1, n, bgs_target,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TLONGLONG, 8, offset + 1, 1, n, mws_target,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TDOUBLE, 9, offset + 1, 1, n, ra, &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TDOUBLE, 10, offset + 1, 1, n, dec, &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TDOUBLE, 11, offset + 1, 1, n, x_focal,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TDOUBLE, 12, offset + 1, 1, n, y_focal,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TSTRING, 13, offset + 1, 1, n, bn_tmp,
                           &status);
            fits_report_error(stderr, status);
            fits_write_col(fptr, TINT, 14, offset + 1, 1, n, fibermask,
                           &status);
            fits_report_error(stderr, status);
        }
        offset += n;
    }
    // PotentialFiberMap table.  We have only one column, so it is safe
    // from a performance perspective to write the whole thing.
    // NO  12/15/16
    strcpy(ttype[0], "POTENTIALTARGETID");
    strcpy(tform[0], "K");
    strcpy(tunit[0], "");
    strcpy(extname, "POTENTIAL_ASSIGNMENTS");
    ret = fits_create_tbl(fptr, BINARY_TBL,
                          potentialtargetid.size(), 1, ttype, tform, tunit, extname, &status);
    fits_report_error(stderr, status);
    fits_write_col(fptr, TLONGLONG, 1, 1, 1, potentialtargetid.size(),
                   &(potentialtargetid[0]), &status);
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);
    for ( size_t c = 0; c < ncols; ++c ) {
        free ( ttype[c] );
        free ( tform[c] );
        free ( tunit[c] );
    }
    free ( ttype );
    free ( tform );
    free ( tunit );
    return;
}
