#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <utils.h>


template <class Ttype> class KDtree {
    private:
        struct TreeNode {
            std::vector <Ttype>  objs;
            double split, weight;
            int axis, parent, left_child, right_child;
        };
        std::vector <TreeNode> tree;
        std::vector <double>   bndbox;
        int Ndim, root, np;
        int partition (std::vector <Ttype> & P, const int left,
                       const int right, const int pivot,
                       const int idim) {
            // This is the partition step of the quickselect algorithm, it is a
            // helper algorithm.
            Ttype tmp;
            double pval = P[pivot].pos[idim];
            tmp = P[P.size() - 1];
            P[P.size() - 1] = P[pivot];
            P[pivot] = tmp;
            int store = left;
            for (int i = left; i <= right; ++i) {
                if (P[i].pos[idim] < pval) {
                    tmp = P[store];
                    P[store] = P[i];
                    P[i] = tmp;
                    store++;
                }
            }
            tmp = P[P.size() - 1];
            P[P.size() - 1] = P[store];
            P[store] = tmp;
            return (store);
        }

        void quickselect (std::vector <Ttype> & P, const int idim,
                          const int mid) {
            // Arranges the array so that objects less than the median are on
            // the
            // left and objects greater than the median are on the right.
            int l = 0;
            int r = P.size() - 1;
            if (l != r) {
                for (;; ) {
                    // rnc 3/28/16 kill randomization
                    // int pivot=(int)(l+drand48()*(r-l+1));
                    int pivot = (int)(l + 0.5 * (r - l + 1) );
                    pivot = partition(P, l, r, pivot, idim);
                    if (pivot == mid) {
                        break;
                    } else if (mid < pivot) {
                        r = pivot - 1;
                    } else {
                        l = pivot + 1;
                    }
                }
            }
        }

        void buildtree (std::vector <Ttype> & P, const int idim) {
            // The main method used to (recursively) build the tree.
            // This uses a post-order "scheme".
            const int Nmax = 32;
            if (P.size() > Nmax) {
                // Will split this node.
                // Select along the relevant axis.
                int median = P.size() / 2;
                quickselect(P, idim % Ndim, median);
                // Make two copies: this is memory expensive, but simple and
                // handles sorts colliding trivially.
                std::vector <Ttype> left, right;
                try {
                    left.assign( P.begin(), P.begin() + median);
                    right.assign(P.begin() + median, P.end() );
                } catch (std::exception & e) {
                    fiberassign::myexception(e);
                }
                // Set up the entry for the split point, this recursively puts
                // in
                // the left and right children--note I am always added last.
                struct TreeNode cur;
                cur.split = P[median].pos[idim];
                cur.axis  = idim;
                if (left.size() > 0) {
                    buildtree(left, (idim + 1) % Ndim);
                    cur.left_child = tree.size() - 1;
                } else {
                    cur.left_child = -1;
                }
                if (right.size() > 0) {
                    buildtree(right, (idim + 1) % Ndim);
                    cur.right_child = tree.size() - 1;
                } else {
                    cur.right_child = -1;
                }
                try {
                    tree.push_back(cur);
                } catch (std::exception & e) {
                    fiberassign::myexception(e);
                }
            } else {
                // Don't split this.
                struct TreeNode cur;
                cur.split = 0;
                cur.axis  = -1;
                try {
                    cur.objs.assign(P.begin(), P.end() );
                } catch (std::exception & e) {
                    fiberassign::myexception(e);
                }
                cur.left_child  = -1;
                cur.right_child = -1;
                try {
                    tree.push_back(cur);
                } catch (std::exception & e) {
                    fiberassign::myexception(e);
                }
            }
        }

        void updatebndry (const int nodeptr) {
            // Uses the parent boundary and current split point to update the
            // boundary.
            int par = tree[nodeptr].parent;
            for (int jdim = 0; jdim < Ndim; ++jdim) {
                bndbox[2 * (Ndim * nodeptr + jdim) + 0] =
                    bndbox[2 * (Ndim * par + jdim) + 0]; // Min
                bndbox[2 * (Ndim * nodeptr + jdim) + 1] =
                    bndbox[2 * (Ndim * par + jdim) + 1]; // Max
            }
            int j = tree[par].axis;
            if (nodeptr == tree[par].left_child) {
                bndbox[2 * (Ndim * nodeptr + j) + 1] = tree[par].split;
            } else {
                bndbox[2 * (Ndim * nodeptr + j) + 0] = tree[par].split;
            }
            if (tree[nodeptr].left_child > -1) {
                updatebndry(tree[nodeptr].left_child);
            }
            if (tree[nodeptr].right_child > -1) {
                updatebndry(tree[nodeptr].right_child);
            }
        }

        void setbndry (const std::vector <double> & xmin,
                       const std::vector <double> & xmax) {
            // Sets the bounding box for each tree node.  This is dynamically
            // allocated as one big double array, starting with the root node,
            // which explains the odd format.
            try {
                bndbox.resize(2 * Ndim * tree.size() );
            } catch (std::exception & e) {
                fiberassign::myexception(e);
            }
            for (int jdim = 0; jdim < Ndim; ++jdim) {
                bndbox[2 * (Ndim * root + jdim) + 0] = xmin[jdim];
                bndbox[2 * (Ndim * root + jdim) + 1] = xmax[jdim];
            }
            if (tree[root].left_child > -1) {
                updatebndry(tree[root].left_child);
            }
            if (tree[root].right_child > -1) {
                updatebndry(tree[root].right_child);
            }
        }

        void setparent (const int nodeptr) {
            if (tree[nodeptr].left_child > -1) {
                tree[tree[nodeptr].left_child].parent = nodeptr;
                setparent(tree[nodeptr].left_child);
            }
            if (tree[nodeptr].right_child > -1) {
                tree[tree[nodeptr].right_child].parent = nodeptr;
                setparent(tree[nodeptr].right_child);
            }
        }

        std::vector <int64_t> near (const int n, const double pos[],
                                const double rmin, const double rmax) const {
            // Returns a vector of object ids near pos within node n.
            std::vector <int64_t> rtn;
            std::vector <Ttype> nil;
            _near_internal(n, pos, rmin, rmax, rtn, nil, false);
            return rtn;
        }

        std::vector <Ttype> near_with_data (const int n, const double pos[],
                                const double rmin, const double rmax) const {
            std::vector <int64_t> nil;
            std::vector <Ttype> rtn;
            _near_internal(n, pos, rmin, rmax, nil, rtn, true);
            return rtn;
        }

        void _near_internal (const int n, const double pos[],
                             const double rmin, const double rmax,
                             std::vector <int64_t> & rtn_id,
                             std::vector <Ttype> & rtn_obj,
                             bool save_obj) const {
            std::vector <int64_t> retlist;
            const double rmax2 = rmax * rmax, rmin2 = rmin * rmin;
            double mindst2 = 0;
            for (int jdim = 0; jdim < Ndim; ++jdim) {
                double mindst = 0;
                if (bound(n, jdim, 0) - pos[jdim] > mindst) {
                    mindst = bound(n, jdim, 0) - pos[jdim];
                }
                if (pos[jdim] - bound(n, jdim, 1) > mindst) {
                    mindst = pos[jdim] - bound(n, jdim, 1);
                }
                if ( ( pos[jdim] >= bound(n, jdim, 0) ) &&
                     ( pos[jdim] <= bound(n, jdim, 1) ) ) {
                    mindst = 0;
                }
                mindst2 += mindst * mindst;
            }
            if (mindst2 < rmax2) {
                double maxdst2 = 0;
                for (int jdim = 0; jdim < Ndim; ++jdim) {
                    double maxdst = 0;
                    if (bound(n, jdim, 1) - pos[jdim] > maxdst) {
                        maxdst = bound(n, jdim, 1) - pos[jdim];
                    }
                    if (pos[jdim] - bound(n, jdim, 0) > maxdst) {
                        maxdst = pos[jdim] - bound(n, jdim, 0);
                    }
                    maxdst2 += maxdst * maxdst;
                }
                if (maxdst2 > rmin2) {
                    if ( ( maxdst2 < rmax2) && ( mindst2 > rmin2) &&
                         ( tree[n].objs.size() > 0) ) {
                        // Contents of cell are entirely within range.
                        size_t orig_size;
                        try {
                            if (!save_obj) {
                                orig_size = rtn_id.size();
                                rtn_id.resize(orig_size + tree[n].objs.size());
                            }
                        } catch (std::exception & e) {
                            fiberassign::myexception(e);
                        }
                        if (save_obj)
                            rtn_obj.insert(rtn_obj.end(),
                                           tree[n].objs.begin(), tree[n].objs.end());
                        else {
                            for (size_t i = 0, iout=orig_size; i < tree[n].objs.size(); ++i, ++iout) {
                                rtn_id[iout] = tree[n].objs[i].id;
                            }
                        }
                        return;
                    } else if (nobj(n) > 0) {
                        // Have lists of objects ... do the brute force
                        // calculation.
                        for (int i = 0; i < nobj(n); ++i) {
                            double r2 = 0;
                            for (int jdim = 0; jdim < Ndim; ++jdim) {
                                r2 += (posn(n, i, jdim) - pos[jdim]) *
                                    (posn(n, i, jdim) - pos[jdim]);
                            }
                            if ( ( r2 < rmax2) && ( r2 > rmin2) ) {
                                try {
                                    if (save_obj)
                                        rtn_obj.push_back(tree[n].objs[i]);
                                    else
                                        rtn_id.push_back(tree[n].objs[i].id);
                                } catch (std::exception & e) {
                                    fiberassign::myexception(e);
                                }
                            }
                        }
                    } else {
                        // Traverse down the hierarchy.
                        int m;
                        if ( (m = get_left_child(n) ) >= 0) {
                            _near_internal(m, pos, rmin, rmax, rtn_id, rtn_obj, save_obj);
                        }
                        if ( (m = get_right_child(n) ) >= 0) {
                            _near_internal(m, pos, rmin, rmax, rtn_id, rtn_obj, save_obj);
                        }
                    }
                }
            }
        }

    public:
        KDtree (const std::vector <Ttype> & Pin, const int Ndim0 = 3) {
            // Create and load the tree -- the elements of P are sorted
            // repeatedly.
            // Ttype should have a member .pos[].
            Ndim = Ndim0;
            // Copy in case want to preserve Pin.
            std::vector <Ttype> P = Pin;
            np = P.size();
            // Build the tree recursively and initialize the parent fields.
            buildtree(P, 0);
            root = tree.size() - 1;
            tree[root].parent = root;
            setparent(root);
            // Set the bounding box(es).
            std::vector <double> xmin, xmax;
            try {
                xmin.resize(Ndim);
                xmax.resize(Ndim);
            } catch (std::exception & e) {
                fiberassign::myexception(e);
            }
            for (int jdim = 0; jdim < Ndim; ++jdim) {
                xmin[jdim] = 1e30;
                xmax[jdim] = -1e30;
            }
            for (int nn = 0; nn < np; ++nn) {
                for (int jdim = 0; jdim < Ndim; ++jdim) {
                    if (xmin[jdim] > P[nn].pos[jdim]) {
                        xmin[jdim] = P[nn].pos[jdim];
                    }
                    if (xmax[jdim] < P[nn].pos[jdim]) {
                        xmax[jdim] = P[nn].pos[jdim];
                    }
                }
            }
            setbndry(xmin, xmax);
        }

        ~KDtree (void) {};  // Do nothing.
        void print_stats () const {
            // Print some stats.
            std::cout << "Tree has Ndim=" << Ndim << ", np=" << np <<
                ", tlen=" << tree.size() << std::endl;
        }

        void print_tree () const {
            // Print the positions of some particles and tree nodes.
            for (int n = 0; n < tree.size(); ++n) {
                if (tree[n].objs.size() > 0) {
                    std::cout << "Tree node " << n << " has " <<
                        tree[n].objs.size() << " objects." << std::endl;
                    for (int i = 0; i < tree[n].objs.size(); ++i) {
                        std::cout << std::fixed << std::setw(10) <<
                            std::setprecision(2) << tree[n].objs[i].pos[0] <<
                            " " << std::fixed << std::setw(10) <<
                            std::setprecision(2) << tree[n].objs[i].pos[1] <<
                            " " << std::fixed << std::setw(10) <<
                            std::setprecision(2) << tree[n].objs[i].pos[2] <<
                            std::endl;
                    }
                } else {
                    std::cout << "Tree node " << n << " has axis=" <<
                        tree[n].axis << " and split " << tree[n].split <<
                        " and left=" << tree[n].left_child << " and right=" <<
                        tree[n].right_child << std::endl;
                }
            }
        }

        int getNdim () const {
            // Returns the dimension of the space.
            return (Ndim);
        }

        int getRoot () const {
            // Returns the root node.
            return (root);
        }

        double bound (const int n, const int idim, const int uplow) const {
            // Returns the lower or upper edge of the idim'th element of the
            // bounding box of node "node".
            return (bndbox[2 * (Ndim * n + idim) + uplow]);
        }

        float posn (const int n, const int m, const int idim) const {
            // Returns the position of the m'th object in the n'th node of the
            // cell.  If idim==Ndim, returns the weight, else the component.
            return (tree[n].objs[m].pos[idim]);
        }

        int nobj (const int n) const {
            // Returns the number of objects in cell n.
            return (tree[n].objs.size() );
        }

        int get_left_child (const int n) const {
            return (tree[n].left_child);
        }

        int get_right_child (const int n) const {
            return (tree[n].right_child);
        }

        std::vector <int64_t> near (const double pos[], const double rmin,
                                const double rmax) const {
            return (near(root, pos, rmin, rmax) );
        }

        std::vector <Ttype> near_with_data (const double pos[], const double rmin,
                                const double rmax) const {
            return (near_with_data(root, pos, rmin, rmax) );
        }

};


#endif
