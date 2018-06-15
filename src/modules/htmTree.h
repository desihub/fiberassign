#ifndef _HTMTREE_H_
#define _HTMTREE_H_

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
/*

   A hierarchical triangular mesh (quad)tree structure for returning
   lists of object ids within a given radius of a certain unit vector.

   Author:      Martin White    (UCB)
   Written: 20-May-2014
   Modified:    21-May-2014
 */
extern void myexit (const int flag);
extern void myexception (const std::exception & e);
template <class Ttype> class htmTree {
    protected:
        struct triangle {
            // vertices and center nhat components.
            double v[9], c[3];
            // an angular size.
            double syze;
        };
        struct treenode {
            struct triangle t;
            std::vector <int> lst;
            int indx, child[4];
        };
        static const int MaxList = 1024;
        double min_node_size;
        std::vector <treenode> T;
        void unitize (double vv[]) {
            // Divides vv by its 3-length to make it a unit vector.
            double ll = 0;
            for (int i = 0; i < 3; i++) {
                ll += vv[i] * vv[i];
            }
            ll = sqrt(ll);
            if (ll <= 0) {
                std::cerr << "Vector length 0!" << std::endl;
                myexit(1);
            }
            for (int i = 0; i < 3; i++) {
                vv[i] /= ll;
            }
        }

        void set_center_size (struct triangle & t) {
            // Given the vertices, sets the center -- just the normalized
            // average --
            // and the triangle "size".
            for (int i = 0; i < 3; ++i) {
                t.c[i] = (t.v[3 * 0 + i] + t.v[3 * 1 + i]
                    + t.v[3 * 2 + i]) / 3.;
            }
            unitize(t.c);
            double cmin = 1;
            for (int i = 0; i < 3; ++i) {
                double dotv = 0;
                for (int j = 0; j < 3; ++j) {
                    dotv += t.c[j] * t.v[3 * i + j];
                }
                if (dotv < cmin) cmin = dotv;
            }
            if (cmin < 1) {
                t.syze = acos(cmin);
            } else {
                t.syze = 0;
            }
        }

        bool in_triangle (const Ttype & P, const struct triangle tt) const {
            // True if P lies in triangle tt.
            if ( P.nhat[0] *
                 (tt.v[3 * 0 + 1] * tt.v[3 * 1 + 2] - tt.v[3 * 0 + 2] *
                  tt.v[3 * 1 + 1]) + P.nhat[1] *
                 (tt.v[3 * 0 + 2] * tt.v[3 * 1 + 0] - tt.v[3 * 0 + 0] *
                  tt.v[3 * 1 + 2]) + P.nhat[2] *
                 (tt.v[3 * 0 + 0] * tt.v[3 * 1 + 1] - tt.v[3 * 0 + 1] *
                  tt.v[3 * 1 + 0]) > 0) return (false);
            if ( P.nhat[0] *
                 (tt.v[3 * 1 + 1] * tt.v[3 * 2 + 2] - tt.v[3 * 1 + 2] *
                  tt.v[3 * 2 + 1]) + P.nhat[1] *
                 (tt.v[3 * 1 + 2] * tt.v[3 * 2 + 0] - tt.v[3 * 1 + 0] *
                  tt.v[3 * 2 + 2]) + P.nhat[2] *
                 (tt.v[3 * 1 + 0] * tt.v[3 * 2 + 1] - tt.v[3 * 1 + 1] *
                  tt.v[3 * 2 + 0]) > 0) return (false);
            if ( P.nhat[0] *
                 (tt.v[3 * 2 + 1] * tt.v[3 * 0 + 2] - tt.v[3 * 2 + 2] *
                  tt.v[3 * 0 + 1]) + P.nhat[1] *
                 (tt.v[3 * 2 + 2] * tt.v[3 * 0 + 0] - tt.v[3 * 2 + 0] *
                  tt.v[3 * 0 + 2]) + P.nhat[2] *
                 (tt.v[3 * 2 + 0] * tt.v[3 * 0 + 1] - tt.v[3 * 2 + 1] *
                  tt.v[3 * 0 + 0]) > 0) return (false); return (true);
        }

        int which_child (const Ttype & P, const int cur) const {
            // Figure out which of the 4 children of "cur" this point lies in.
            int c = -1;
            for (int ichild = 0; c == -1 && ichild < 4; ichild++) {
                if (in_triangle(P, T[T[cur].child[ichild]].t) ) c = ichild;
            }
            if (c == -1) {
                // Some kind of round-off error must have
                // occurred for me to get here.  Use brute force.
                double cmax = -1;
                for (int ichild = 0; ichild < 4; ++ichild) {
                    double cdot = 0;
                    for (int i = 0; i < 3; ++i) {
                        cdot += P.nhat[i] * T[T[cur].child[ichild]].t.c[i];
                    }
                    if (cdot > cmax) {
                        cmax = cdot;
                        c    = ichild;
                    }
                }
            }
            return (c);
        }

        void create_children (const int nn) {
            // Creates the 4 children of node nn, having them point to -2.
            struct treenode tn;
            // Child 0.
            // Vertex 0 is the same.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 0 + i] = T[nn].t.v[3 * 0 + i];
            }
            unitize(&tn.t.v[0]);
            // Vertex 1 is midway between 0 & 1.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 1 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 1 + i]) / 2.;
            }
            unitize(&tn.t.v[3]);
            // Vertex 2 is midway between 0 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 2 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[6]);
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T[nn].child[0] = T.size();
            try {
                T.push_back(tn);
            } catch (std::exception & e) {
                myexception(e);
            }
            // Child 1.
            // Vertex 0 is midway between 0 & 1.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 0 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 1 + i]) / 2.;
            }
            unitize(&tn.t.v[0]);
            // Vertex 1 is the same
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 1 + i] = T[nn].t.v[3 * 1 + i];
            }
            unitize(&tn.t.v[3]);
            // Vertex 2 is midway between 1 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 2 + i] =
                    (T[nn].t.v[3 * 1 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[6]);
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T[nn].child[1] = T.size();
            try {
                T.push_back(tn);
            } catch (std::exception & e) {
                myexception(e);
            }
            // Child 2.
            // Vertex 0 is midway between 0 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 0 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[0]);
            // Vertex 1 is midway between 1 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 1 + i] =
                    (T[nn].t.v[3 * 1 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[3]);
            // Vertex 2 is the same.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 2 + i] = T[nn].t.v[3 * 2 + i];
            }
            unitize(&tn.t.v[6]);
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T[nn].child[2] = T.size();
            try {
                T.push_back(tn);
            } catch (std::exception & e) {
                myexception(e);
            }
            // Child 3.
            // Vertex 0 is midway between 1 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 0 + i] =
                    (T[nn].t.v[3 * 1 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[0]);
            // Vertex 1 is midway between 0 & 2.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 1 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 2 + i]) / 2.;
            }
            unitize(&tn.t.v[3]);
            // Vertex 2 is midway between 0 & 1.
            for (int i = 0; i < 3; ++i) {
                tn.t.v[3 * 2 + i] =
                    (T[nn].t.v[3 * 0 + i] + T[nn].t.v[3 * 1 + i]) / 2.;
            }
            unitize(&tn.t.v[6]);
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T[nn].child[3] = T.size();
            try {
                T.push_back(tn);
            } catch (std::exception & e) {
                myexception(e);
            }
        }

        void set_base (void) {
            // First set up the basic 8 triangles on the sphere and have
            // them point to nothing.  This also clears the tree.
            T.clear();
            struct treenode tn;
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] =  1;
            tn.t.v[3] =  0;
            tn.t.v[4] =  1;
            tn.t.v[5] =  0;
            tn.t.v[6] =  1;
            tn.t.v[7] =  0;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] =  1;
            tn.t.v[3] = -1;
            tn.t.v[4] =  0;
            tn.t.v[5] =  0;
            tn.t.v[6] =  0;
            tn.t.v[7] =  1;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] =  1;
            tn.t.v[3] =  0;
            tn.t.v[4] = -1;
            tn.t.v[5] =  0;
            tn.t.v[6] = -1;
            tn.t.v[7] =  0;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] =  1;
            tn.t.v[3] =  1;
            tn.t.v[4] =  0;
            tn.t.v[5] =  0;
            tn.t.v[6] =  0;
            tn.t.v[7] = -1;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] = -1;
            tn.t.v[3] =  1;
            tn.t.v[4] =  0;
            tn.t.v[5] =  0;
            tn.t.v[6] =  0;
            tn.t.v[7] =  1;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] = -1;
            tn.t.v[3] =  0;
            tn.t.v[4] =  1;
            tn.t.v[5] =  0;
            tn.t.v[6] = -1;
            tn.t.v[7] =  0;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] = -1;
            tn.t.v[3] = -1;
            tn.t.v[4] =  0;
            tn.t.v[5] =  0;
            tn.t.v[6] =  0;
            tn.t.v[7] = -1;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
            tn.t.v[0] =  0;
            tn.t.v[1] =  0;
            tn.t.v[2] = -1;
            tn.t.v[3] =  0;
            tn.t.v[4] = -1;
            tn.t.v[5] =  0;
            tn.t.v[6] =  1;
            tn.t.v[7] =  0;
            tn.t.v[8] =  0;
            set_center_size(tn.t);
            for (int i = 0; i < 4; ++i) {
                tn.child[i] = -1;
            }
            tn.indx = -2;
            T.push_back(tn);
        }

        void add_object (const int i, const std::vector <Ttype> & P) {
            // Adds on object to the tree.
            bool notput = true;
            for (int ibase = 0; ibase < 8 && notput; ++ibase) {
                if (in_triangle(P[i], T[ibase].t) ) {
                    int cur = ibase;
                    while (notput) {
                        if (T[cur].indx == -2) {
                            T[cur].indx = i;
                            notput = false;
                        } else if ( ( T[cur].indx == -1) &&
                                    ( T[cur].lst.size() > 0) ) {
                            // I am an internal node with a list, simply
                            // append.
                            try {
                                T[cur].lst.push_back(i);
                            } catch (std::exception & e) {
                                myexception(e);
                            }
                            notput = false;
                        } else if ( ( T[cur].indx == -1) &&
                                    ( T[cur].lst.size() == 0) ) {
                            // Current node is an "internal" node, it has
                            // children
                            // but nothing in its list.
                            int c = which_child(P[i], cur);
                            cur  = T[cur].child[c];
                        } else if (T[cur].indx >= 0) {
                            // Current node points to a particle--need to
                            // refine.
                            int pid = T[cur].indx;
                            if (T[cur].t.syze > min_node_size) {
                                // Should subdivide.
                                create_children(cur);
                                T[cur].indx = -1;
                                int c = which_child(P[pid], cur);
                                T[T[cur].child[c]].indx = pid;
                                c = which_child(P[i], cur);
                                cur = T[cur].child[c];
                            } else {
                                // Promote to a terminal node and put in
                                // particle.
                                T[cur].indx = -1;
                                try {
                                    T[cur].lst.push_back(pid);
                                    T[cur].lst.push_back(i);
                                } catch (std::exception & e) {
                                    myexception(e);
                                }
                                notput = false;
                            }
                        }
                    }
                }
            }
        }

    public:
        ~htmTree () {} // Do nothing.

        void load (const std::vector <Ttype> & P) {
            // Load an existing tree with objects.
            set_base();
            for (int i = 0; i < P.size(); ++i) {
                add_object(i, P);
            }
        }

        htmTree () {
            // Set up an empty tree.
            min_node_size = 1e-3;
            set_base();
        }

        htmTree (std::vector <Ttype> & P, const double minsize = 0.001) {
            // Create a tree and load it with objects P.
            min_node_size = minsize;
            load(P);
        }

        void stats (void) const {
            // Prints some (possibly) useful statistics.
            double minsize = 1e10, smlsize = 100.;
            int maxlen =  -1, totlen = 0;
            for (int i = 0; i < T.size(); i++) {
                double lsiz = T[i].t.syze;
                // Need to do this for some reason.
                int llen = T[i].lst.size();
                if (minsize > lsiz) {
                    minsize = lsiz;
                }
                if (maxlen  < llen) {
                    maxlen = llen;
                    smlsize = lsiz;
                }
                totlen += llen;
            }
            std::cout << "# Tree has " << T.size() << " nodes." << std::endl;
            std::cout << "# Memory load is " << T.size() *
                sizeof(struct treenode) / 1024. / 1024. <<
                " MB (excluding lists)." << std::endl;
            std::cout << "# Longest list length is " << maxlen << std::endl;
            std::cout << "# Occurs for node of size " << smlsize << std::endl;
            std::cout << "# Total cumulative list length is " << totlen <<
                std::endl;
            std::cout << "# Minimum node size is " << std::scientific <<
                minsize << " (min_node_size=" << min_node_size << ")." <<
                std::endl;
        }

        void which_base (const Ttype & P) {
            // Prints which base(s) point P lies in.
            for (int ibase = 0; ibase < 8; ++ibase) {
                if (in_triangle(P, T[ibase].t) ) {
                    std::cout << std::fixed << std::setw(8) <<
                        std::setprecision(4) << P.nhat[0] << std::fixed <<
                        std::setw(8) << std::setprecision(4) << P.nhat[1] <<
                        std::fixed << std::setw(8) << std::setprecision(4) <<
                        P.nhat[2] << std::fixed << std::setw(8) << ibase <<
                        std::endl;
                }
            }
        }

        std::vector <int> near (const std::vector <Ttype> & P,
                                const double nhat[], const double rad) const {
            // Returns a list (vector) of indices of objects "near" nhat.
            // Distances refer to \theta, not 2sin(theta/2).
            const double cosr = cos(rad);
            std::vector <int> clist, nbr;
            for (int i = 0; i < 8; ++i) {
                clist.push_back(i);
            }
            while (clist.size() > 0) {
                int cur = clist.back();
                clist.pop_back();
                if (T[cur].indx >= 0) {
                    // I am a leaf pointing at a particle.
                    int indx = T[cur].indx;
                    double dotv = 0, dst = 0;
                    for (int i = 0; i < 3; i++) {
                        dotv += nhat[i] * P[indx].nhat[i];
                    }
                    if (dotv < 1) {
                        dst = acos(dotv);
                    }
                    if (dst < rad) {
                        try {
                            nbr.push_back(indx);
                        } catch (std::exception & e) {
                            myexception(e);
                        }
                    }
                } else if ( ( T[cur].indx == -1) && ( T[cur].lst.size() > 0) ) {
                    // I am a leaf node with a list.
                    for (int j = 0; j < T[cur].lst.size(); ++j) {
                        double dotv = 0, dst = 0;
                        for (int i = 0; i < 3; i++) {
                            dotv += nhat[i] * P[T[cur].lst[j]].nhat[i];
                        }
                        if (dotv < 1) {
                            dst = acos(dotv);
                        }
                        if (dst < rad) {
                            try {
                                nbr.push_back(T[cur].lst[j]);
                            } catch (std::exception & e) {
                                myexception(e);
                            }
                        }
                    }
                } else if ( ( T[cur].indx == -1) &&
                            ( T[cur].lst.size() == 0) ) {
                    // I am an internal node, with children.
                    // Compute "distance" from nhat to center of triangle.
                    double dotv = 0;
                    for (int i = 0; i < 3; i++) {
                        dotv += nhat[i] * T[cur].t.c[i];
                    }
                    // and subtract "size" of triangle.
                    // If the smallest distance lies within search, open ...
                    double dmin = acos(dotv) - T[cur].t.syze;
                    if ( ( dotv >= 1) || ( dmin < rad) ) {
                        for (int i = 0; i < 4; ++i) {
                            clist.push_back(T[cur].child[i]);
                        }
                    }
                }
            }
            return (nbr);
        }
};


#endif
