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
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <sys/time.h>
#include <map>
#include <cstring>
#include <cstdint>
#include "misc.h"
#include "feat.h"
#include "structs.h"
#include "collision.h"
#include "fitsio.h"

// targets
// -----------------------------------------------------------------------
MTL read_MTLfile (str readfile, const Feat & F, long SS, long SF) {
    // reads fits files, specifically mtl, standard stars, sky fibers
    str s = readfile;
    MTL M;
    std::string buf;
    const char * fname;
    fname = s.c_str();
    std::ifstream fs(fname);
    int ii;
    fitsfile * fptr;
    int status = 0, anynulls;
    int hdutype;
    int nkeys;
    int hdupos;
    long nrows;
    long nkeep;
    int ncols;
    long long * targetid;
    long * desi_target;
    long * bgs_target;
    long * mws_target;
    long starmask = SS;
    int * numobs;
    int * priority;
    double * ra;
    double * dec;
    int colnum;
    char * * brickname;
    uint16_t * obsconditions;
    double * subpriority;
    fprintf(stdout, "star mask %ld\n", starmask);
    // General purpose output stream for exceptions
    std::ostringstream o;
    // Check that input file exists and is readable by cfitsio
    std::cout << "Finding file: " << fname << std::endl;
    int file_exists;
    fits_file_exists(fname, &file_exists, &status);
    std::ostringstream exists_str;
    exists_str << "(CFITSIO file_exists code: " << file_exists << ")";
    // Throw exceptions for failed read, see cfitsio docs
    if (!file_exists) {
        switch (file_exists) {
            case -1: o << "Input MTL file must be a disk file: " << fname <<
                " " << exists_str.str();
                throw std::runtime_error(o.str().c_str() );
            case  0: o << "Could not find MTL input file: " << fname << " " <<
                exists_str.str();
                throw std::runtime_error(o.str().c_str() );
            case  2: o << "Cannot handle zipped MTL input file: " << fname <<
                " " << exists_str.str();
                throw std::runtime_error(o.str().c_str() );
        }
    }
    std::cout << "Found MTL input file: " << fname << std::endl;
    if (!fits_open_file(&fptr, fname, READONLY, &status) ) {
        std::cout << "Reading MTL input file " << fname << std::endl;
        if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) myexit(status);
        fits_get_hdrspace(fptr, &nkeys, NULL, &status);
        fits_get_hdu_num(fptr, &hdupos);
        // Get the HDU type
        fits_get_hdu_type(fptr, &hdutype, &status);
        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);
        fflush(stdout);
        if (!(targetid = (long long *)malloc(nrows * sizeof(long long) ) ) ) {
            fprintf(stderr, "problem with targetid allocation\n");
            myexit(1);
        }
        if (!(desi_target = (long *)malloc(nrows * sizeof(long) ) ) ) {
            fprintf(stderr, "problem with desi_target allocation\n");
            myexit(1);
        }
        if (!(mws_target = (long *)malloc(nrows * sizeof(long) ) ) ) {
            fprintf(stderr, "problem with mws_target allocation\n");
            myexit(1);
        }
        if (!(bgs_target = (long *)malloc(nrows * sizeof(long) ) ) ) {
            fprintf(stderr, "problem with bgs_target allocation\n");
            myexit(1);
        }
        if (!(numobs = (int *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with numobs allocation\n");
            myexit(1);
        }
        if (!(priority = (int *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with priority allocation\n");
            myexit(1);
        }
        if (!(subpriority = (double *)malloc(nrows * sizeof(double) ) ) ) {
            fprintf(stderr, "problem with priority allocation\n");
            myexit(1);
        }
        if (!(ra = (double *)malloc(nrows * sizeof(double) ) ) ) {
            fprintf(stderr, "problem with ra allocation\n");
            myexit(1);
        }
        if (!(dec = (double *)malloc(nrows * sizeof(double) ) ) ) {
            fprintf(stderr, "problem with dec allocation\n");
            myexit(1);
        }
        if (!(obsconditions = (uint16_t *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with obsconditions allocation\n");
            myexit(1);
        }
        if (!(brickname = (char **)malloc(nrows * sizeof(char *) ) ) ) {
            fprintf(stderr, "problem with brickname allocation\n");
            myexit(1);
        }
        for (ii = 0; ii < nrows; ii++) {
            if (!(brickname[ii] = (char *)malloc(9 * sizeof(char) ) ) ) {
                fprintf(stderr, "problem with brickname allocation\n");
                myexit(1);
            }
        }

        // ----- TARGETID
        // find which column contains the TARGETID values
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"TARGETID", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding TARGETID column\n");
            myexit(status);
        }
        long frow, felem, nullval;
        frow = 1;
        felem = 1;
        nullval = -99.;
        if (fits_read_col(fptr, TLONGLONG, colnum, frow, felem, nrows,
                          &nullval, targetid, &anynulls, &status) ) {
            fprintf(stderr, "error reading TARGETID column\n");
            myexit(status);
        }

        // ----- RA
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"RA", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding RA column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nrows, &nullval,
                          ra, &anynulls, &status) ) {
            fprintf(stderr, "error reading RA column\n");
            myexit(status);
        }

        // ----- DEC
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"DEC", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding DEC column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nrows, &nullval,
                          dec, &anynulls, &status) ) {
            fprintf(stderr, "error reading DEC column\n");
            myexit(status);
        }

        // ----- Target mask bits
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"DESI_TARGET", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding DESI_TARGET column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TLONG, colnum, frow, felem, nrows, &nullval,
                          desi_target, &anynulls, &status) ) {
            fprintf(stderr, "error reading DESI_TARGET column\n");
            myexit(status);
        }
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"MWS_TARGET", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding MWS_TARGET column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TLONG, colnum, frow, felem, nrows, &nullval,
                          mws_target, &anynulls, &status) ) {
            fprintf(stderr, "error reading MWS_TARGET column\n");
            myexit(status);
        }
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"BGS_TARGET", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding BGS_TARGET column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TLONG, colnum, frow, felem, nrows, &nullval,
                          bgs_target, &anynulls, &status) ) {
            fprintf(stderr, "error reading BGS_TARGET column\n");
            myexit(status);
        }

        // OBSCONDITIONS
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"OBSCONDITIONS", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding OBSCONDITIONS column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, USHORT_IMG, colnum, frow, felem, nrows,
                          &nullval, obsconditions, &anynulls, &status) ) {
            fprintf(stderr, "error reading OBSCONDITIONS column\n");
            myexit(status);
        }

        // ----- BRICKNAME
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"BRICKNAME", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding BRICKNAME column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TSTRING, colnum, frow, felem, nrows, &nullval,
                          brickname, &anynulls, &status) ) {
            fprintf(stderr, "error reading BRICKNAME column\n");
            myexit(status);
        }

        // ----- SUBPRIORITY
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"SUBPRIORITY", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding SUBPRIORITY column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nrows, &nullval,
                          subpriority, &anynulls, &status) ) {
            fprintf(stderr, "error reading SUBPRIORITY column\n");
            myexit(status);
        }

        // ----- NUMOBS_MORE
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"NUMOBS_MORE", &colnum,
                             &status) ) {
            // fprintf(stderr, "error finding NUMOBS_MORE column\n");
            // myexit(status);
            std::cout << "NUMOBS_MORE not found ... setting to 1" << std::endl;
            for (int i = 0; i < nrows; i++) {
                numobs[i] = 1;
            }
        } else if (fits_read_col(fptr, TINT, colnum, frow, felem, nrows,
                                 &nullval, numobs, &anynulls, &status) ) {
            fprintf(stderr, "error reading NUMOBS_MORE column\n");
            myexit(status);
        }
        // ----- PRIORITY
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"PRIORITY", &colnum,
                             &status) ) {
            // fprintf(stderr, "error finding PRIORITY column\n");
            // myexit(status);
            std::cout << "PRIORITY not found ... setting to 1" << std::endl;
            for (int i = 0; i < nrows; i++) {
                priority[i] = 1;
            }
        } else if (fits_read_col(fptr, TINT, colnum, frow, felem, nrows,
                                 &nullval, priority, &anynulls, &status) ) {
            fprintf(stderr, "error reading PRIORITY column\n");
            myexit(status);
        }
        nkeep = 0;
        for (ii = 0; ii < nrows; ii++) {
            if ( ( (SS != 0) && ( (desi_target[ii] & starmask) != 0) ) ||
                 (SS == 0) ) {
                nkeep++;
            }
        }
        // count how many rows we will keep and reserve that amount
        // nkeep = nrows;
        printf("Keeping %ld targets within ra/dec ranges\n", nkeep);
        try {
            M.reserve(nkeep);
        } catch (std::exception & e) {
            myexception(e);
        }
        for (ii = 0; ii < nrows; ii++) {
            str xname;
            // make sure ra is between 0 and 360
            if (ra[ii] <   0.) ra[ii] += 360.;
            if (ra[ii] >= 360.) ra[ii] -= 360.;
            if ( ( dec[ii] <= -90.) || ( dec[ii] >= 90.) ) {
                std::cout << "DEC=" << dec[ii] << " out of range reading " <<
                fname << std::endl;
                myexit(1);
            }
            double theta = (90.0 - dec[ii]) * M_PI / 180.;
            double phi   = (ra[ii]        ) * M_PI / 180.;
            struct target Q;
            Q.nhat[0]    = cos(phi) * sin(theta);
            Q.nhat[1]    = sin(phi) * sin(theta);
            Q.nhat[2]    = cos(theta);
            Q.obsconditions = obsconditions[ii];
            // priority not present for sky fibers or standard stars
            Q.t_priority = priority[ii];
            Q.subpriority = subpriority[ii];
            Q.nobs_remain = numobs[ii];
            // need to keep track of this, too
            Q.nobs_done = 0;
            // changed only in update_plan
            Q.once_obs = 0;
            Q.ra = ra[ii];
            Q.dec = dec[ii];
            Q.id = targetid[ii];
            Q.desi_target = desi_target[ii];
            Q.mws_target = mws_target[ii];
            Q.bgs_target = bgs_target[ii];

	    if (SS != 0) {
	      Q.SS = 1;
	    } else {
	      Q.SS = SS;
	    }
            Q.SF = SF;

            // These variables were not initialized elsewhere
            Q.lastpass = 0;
            Q.priority_class = 0;
            strncpy(Q.brickname, brickname[ii], 9);
            if ( ( (SS != -1) && ( (desi_target[ii] & starmask) != 0) ) ||
                 (SS == 0) ) {
                try {
                    M.push_back(Q);
                } catch (std::exception & e) {
                    myexception(e);
                }
                bool in = false;
                for (size_t j = 0; j < M.priority_list.size(); ++j) {
                    if (Q.t_priority == M.priority_list[j]) {
                        in = true;
                    }
                }
                if (!in) {
                    M.priority_list.push_back(Q.t_priority);
                }
            }
        }  // end ii loop over targets

        // Free memory
        free(targetid);
        free(desi_target);
        free(mws_target);
        free(bgs_target);
        free(numobs);
        free(priority);
        free(subpriority);
        free(ra);
        free(dec);
        free(obsconditions);
        for (ii = 0; ii < nrows; ii++) {
            free(brickname[ii]);
        }
        free(brickname);

        // Close file
        fits_close_file(fptr, &status);

        std::sort(M.priority_list.begin(), M.priority_list.end() );
        return (M);
    } else {
        std::ostringstream open_status_str;
        open_status_str << "(CFITSIO open_file status: " << status << ")";
        o << "Problem opening input MTL fits file: " << fname << " " <<
        open_status_str.str();
        throw std::runtime_error(o.str().c_str() );
    }
}

void assign_priority_class (MTL & M) {
    // assign each target to a priority class
    // this needs to be updated
    for (size_t i = 0; i < M.size(); ++i) {
        if (!M[i].SS && !M[i].SF) {
            for (size_t j = 0; j < M.priority_list.size(); ++j) {
                if (M[i].t_priority == M.priority_list[j]) {
                    M[i].priority_class = j;
                }
            }
        }
    }
}

// FP
// ----------------------------------------------------------------------------
// Read the positions of the fibers on each plate.
// need also to get the petal, i.e. spectrometer  rnc 1/16/15  added S
FP  read_fiber_positions (const Feat & F) {
    std::string buf;
    std::ifstream fs(F.fibFile.c_str() );
    if (!fs) {
        // An error occurred opening the file.
        std::cerr << "Unable to open file " << F.fibFile << std::endl;
        myexit(1);
    }
    getline(fs, buf);
    while (fs.eof() == 0 &&
           ( (buf[0] == '#') || (buf.size() == 0) ) ) {
        getline(fs, buf);
    }
    int i(0);
    int petals_pac[] = {0, 1, 2, 7, 8, 9};
    List petals_pacL = initList(petals_pac, 6);
    List inv = inverse(petals_pacL);
    // collection of fibersfound(spectro,petals_pacL))
    FP FibPos;
    fpos fiber_pos;
    printf("before reading positioners \n");
    std::cout.flush();
    while (fs.eof() == 0) {
        double x, y;
        int fiber, location, spectro;
        std::istringstream(buf) >> fiber >> location >> spectro >> x >> y;
        try {
            fiber_pos.fib_num = fiber;
            fiber_pos.location = location;
            fiber_pos.fp_x = x;
            fiber_pos.fp_y = y;
            //int sp =  spectro;
            fiber_pos.spectrom = spectro;
            fiber_pos.coords = dpair(x, y);
        } catch (std::exception & e) {
            myexception(e);
        }
        FibPos.push_back(fiber_pos);
        getline(fs, buf);
        i++;
    }
    fs.close();
    printf("read the positioner file\n");
    int fiber_size = FibPos.size();
    // sort by fiber number
    std::vector <std::pair <int, int> > pairs;
    for (size_t f = 0; f < FibPos.size(); ++f) {
        std::pair <int, int> this_pair (FibPos[f].fib_num, f);
        pairs.push_back(this_pair);
    }
    std::sort(pairs.begin(), pairs.end(), int_pairCompare);
    std::vector <fpos> out;
    for (size_t f = 0; f < FibPos.size(); ++f) {
        out.push_back(FibPos[pairs[f].second]);
    }
    copy(out.begin(), out.end(), FibPos.begin() );
    printf(" sorted by fiber number\n");
    for (int i = 0; i < 10; ++i) {
        printf(" i %d FibPos[i].fib_num %d \n", i, FibPos[i].fib_num);
    }
    // create fibers_of_sp
    FibPos.fibers_of_sp.resize(F.Npetal);
    for (int k = 0; k < fiber_size; k++) {
        FibPos.fibers_of_sp[FibPos[k].spectrom].push_back(k);
    }
    // create table of Neighbors
    for (int i = 0; i < fiber_size; i++) {
        for (int j = 0; j < fiber_size; j++) {
            if (i != j) {
                if (sq(FibPos[i].fp_x - FibPos[j].fp_x) +
                       sq(FibPos[i].fp_y - FibPos[j].fp_y) <
                       sq(F.NeighborRad) ) {
                    FibPos[i].N.push_back(j);
                }
            }
        }
    }
    printf(" made neighbors \n");
    return (FibPos);
}

// FP::FP() {};
int A_less_than_B (int year_A, int month_A, int day_A, int year_B, int month_B,
                   int day_B) {
    // fprintf(stdout, "%d %d %d %d %d %d\n",year_A, month_A, day_A, year_B,
    // month_B, day_B);
    if ( ( (year_A + month_A / 12.0) < (year_B + month_B / 12.0) ) ) return 1;
    if ( ( (year_A + month_A / 12.0) > (year_B + month_B / 12.0) ) ) return 0;
    if ( ( (year_A + month_A / 12.0) == (year_B + month_B / 12.0) ) ) {
        if (day_A < day_B) {
            return 1;
        } else {
            return 0;
        }
    }
    return -1;
}

void read_fiber_status (FP & FibPos, const Feat & F) {
    std::string buf;
    //char date_now[512];
    char date_init[512];
    char date_end[512];
    std::ifstream fs(F.fibstatusFile.c_str() );
    int fiber_size = FibPos.size();
    std::tm input_time = {};
    std::tm init_time = {};
    std::tm end_time = {};
    std::tm current_time = {};
    // current time
    auto now = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now() );
    current_time = *(std::gmtime(&now) );
    // input time
    std::stringstream ss(F.runDate);
    ss >> std::get_time(&input_time, "%Y-%m-%d");
    // if the input time was set in the input file, then override the current
    // machine time with the input time
    if (A_less_than_B(input_time.tm_year, input_time.tm_mon,
                      input_time.tm_mday, 1, 1, 1) ) {
        std::cout << "INPUT time is not set";
    } else {
        current_time = input_time;
    }
    std::cout << "Input Time" <<  std::put_time(&input_time, "%c") << "\n";
    std::cout << "Current Time" <<  std::put_time(&current_time, "%c") << "\n";
    if (!fs) {
        // An error occurred opening the file.
        std::cerr << "Unable to open file " << F.fibstatusFile << std::endl;
        myexit(1);
    }
    getline(fs, buf);
    while (fs.eof() == 0 && ( (buf[0] == '#') || (buf.size() == 0) ) ) {
        getline(fs, buf);
    }
    std::cout.flush();
    int i(0);
    printf("before reading status \n");
    std::cout.flush();
    while (fs.eof() == 0) {
        double x, y;
        int fiber, location, broken, stuck;
        getline(fs, buf);
        std::istringstream(buf) >> fiber >> location >> x >> y >> broken >>
        stuck >> date_init >> date_end;
        fprintf(stdout,
                "Read from fiber status: Fiber_pos %d Location %d Broken %d Stuck %d dates %s %s \n", fiber, location, broken, stuck, date_init,
                date_end);
        std::stringstream si(date_init);
        si >> std::get_time(&init_time, "%Y-%m-%dT%H:%M:%S");
        std::cout << "Init Time for Fiber" <<
        std::put_time(&init_time, "%c") << "\n";
        // std::stringstream ss(date_end);
        std::stringstream se(date_end);
        se >> std::get_time(&end_time, "%Y-%m-%dT%H:%M:%S");
        std::cout << "End Time for Fiber" <<
        std::put_time(&end_time, "%c") << "\n";
        i++;
        // if the current time is within the time interval for the stuck/broken
        // fiber, then make the change
        if (A_less_than_B(init_time.tm_year, init_time.tm_mon,
                          init_time.tm_mday, current_time.tm_year,
                          current_time.tm_mon,
                          current_time.tm_mday) &&
            A_less_than_B(current_time.tm_year, current_time.tm_mon,
                          current_time.tm_mday, end_time.tm_year,
                          end_time.tm_mon,
                          end_time.tm_mday) ) {
            for (int j = 0; j < fiber_size; j++) {
                if (fiber == FibPos[j].fib_num) {
                    if (location == FibPos[j].location) {
                        fprintf(stdout,
                                "Changing fiberastatus entry: Fiber %d Location %d\n", fiber,
                                location);
                        if (broken == 1) {
                            fprintf(stdout, "BROKEN\n");
                            FibPos[j].broken = broken;
                        }
                        if (stuck == 1) {
                            fprintf(stdout, "STUCK\n");
                            FibPos[j].stuck = stuck;
                            FibPos[j].fp_x = x;
                            FibPos[j].fp_y = y;
                        }
                    } else {
                        std::cerr << "Fiber ID Matches But Not Location ID " <<
                        F.fibFile << std::endl;
                        myexit(1);
                    }
                }
            }
        }
    }
    fs.close();
    printf("read status file\n");
}

// plate
// ---------------------------------------------------------------------------
List plate::av_gals_plate (const Feat & F, const MTL & M,
                           const FP & pp) const {
    // list of galaxies available to plate no repetitions
    List gals = initList(F.Ngal);
    List L = initList(0);
    for (int k = 0; k < F.Nfiber; k++) {
        for (size_t i = 0; i < av_gals[k].size(); i++) {
            if (gals[av_gals[k][i]] == 0) {
                gals[av_gals[k][i]] = 1;
                L.push_back(i);
            }
        }
    }
    return L;
}

Plates read_plate_centers (const Feat & F) {
    Plates P, PQ;
    const char * fname;
    /*Variables used to read fits file*/
    fitsfile * fptr;
    int status = 0, anynulls;
    int hdutype;
    int nkeys;
    int hdupos;
    long nrows;
    //long nkeep;
    int ncols;
    int ii;
    uint16_t * obsconditions;
    int * in_desi;
    int * tile_id;
    int tileid;
    int * ipass;
    double * ra;
    double * dec;
    int colnum;
    long frow, felem, nullval;
    frow = 1;
    felem = 1;
    nullval = -99.;
    Time t, time;  // t for global, time for local
    init_time(t);
    // read the strategy file
    // survey_list is list of tiles specified by tileid (arbitrary int) in
    // order of survey
    std::ifstream fsurvey(F.surveyFile.c_str() );
    if (!fsurvey) {
        // An error occurred opening the file.
        std::cerr << "Unable to open file " << F.surveyFile << std::endl;
        myexit(1);
    }
    int survey_tile;
    std::vector <int> survey_list;
    std::string buf;
    while (getline(fsurvey, buf) ) {
        std::istringstream ss(buf);
        if (!(ss >> survey_tile) ) break;
        survey_list.push_back(survey_tile);
    }
    printf(" number of tiles %lu \n", survey_list.size() );
    // NEW
    // read list of tile centers
    // Check that input file exists and is readable by cfitsio
    fname = F.tileFile.c_str();
    std::cout << "Finding file: " << fname << std::endl;
    int file_exists;
    fits_file_exists(fname, &file_exists, &status);
    std::ostringstream exists_str;
    exists_str << "(CFITSIO file_exists code: " << file_exists << ")";
    std::ostringstream o;
    // Throw exceptions for failed read, see cfitsio docs
    if (!file_exists) {
        switch (file_exists) {
            case -1: o << "Input tile centers file must be a disk file: " <<
                fname << " " << exists_str.str();
                throw std::runtime_error(o.str().c_str() );
            case  0: o << "Could not find input tile centers file: " <<
                fname << " " << exists_str.str();
                throw std::runtime_error(o.str().c_str() );
            case  2: o << "Cannot handle zipped tile centers input file: " <<
                fname << " " << exists_str.str();
                throw std::runtime_error(o.str().c_str() );
        }
    }
    std::cout << "Found input tile centers file: " << fname << std::endl;
    if (!fits_open_file(&fptr, fname, READONLY, &status) ) {
        std::cout << "Reading input tile centers file " << fname << std::endl;
        if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) myexit(status);
        fits_get_hdrspace(fptr, &nkeys, NULL, &status);
        fits_get_hdu_num(fptr, &hdupos);
        fits_get_hdu_type(fptr, &hdutype, &status);  /* Get the HDU type */
        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);
        /*
           std::cout << ncols << " columns " << nrows << "nrows" << std::endl;
           std::cout << "HDU " << hdupos << std::endl;
           if (hdutype == ASCII_TBL){
           std::cout << "ASCII TABLE: " << std::endl;
           }else{
           std::cout << "BINARY TABLE: " << std::endl;
           }
         */
        if (!(obsconditions = (uint16_t *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with priority allocation\n");
            myexit(1);
        }
        if (!(ipass = (int *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with ipass allocation\n");
            myexit(1);
        }
        if (!(in_desi = (int *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with priority allocation\n");
            myexit(1);
        }
        if (!(tile_id = (int *)malloc(nrows * sizeof(int) ) ) ) {
            fprintf(stderr, "problem with priority allocation\n");
            myexit(1);
        }
        if (!(ra = (double *)malloc(nrows * sizeof(double) ) ) ) {
            fprintf(stderr, "problem with ra allocation\n");
            myexit(1);
        }
        if (!(dec = (double *)malloc(nrows * sizeof(double) ) ) ) {
            fprintf(stderr, "problem with dec allocation\n");
            myexit(1);
        }

        // ----- RA
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"RA", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding RA column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nrows, &nullval,
                          ra, &anynulls, &status) ) {
            fprintf(stderr, "error reading RA column\n");
            myexit(status);
        }

        // ----- DEC
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"DEC", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding DEC column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nrows, &nullval,
                          dec, &anynulls, &status) ) {
            fprintf(stderr, "error reading DEC column\n");
            myexit(status);
        }

        // ----- IN_DESI
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"IN_DESI", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding IN_DESI column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TINT, colnum, frow, felem, nrows, &nullval,
                          in_desi, &anynulls, &status) ) {
            fprintf(stderr, "error reading IN_DESI column\n");
            myexit(status);
        }

        // ----- OBSCONDITIONS
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"OBSCONDITIONS", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding OBSCONDITIONS column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, USHORT_IMG, colnum, frow, felem, nrows,
                          &nullval, obsconditions, &anynulls, &status) ) {
            fprintf(stderr, "error reading OBSCONDITIONS column\n");
            myexit(status);
        }

        // ----- TILEID
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"TILEID", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding OBSCONDITIONS column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TINT, colnum, frow, felem, nrows, &nullval,
                          tile_id, &anynulls, &status) ) {
            fprintf(stderr, "error reading TILEID column\n");
            myexit(status);
        }

        // ----- PASS
        if ( fits_get_colnum(fptr, CASEINSEN, (char *)"PASS", &colnum,
                             &status) ) {
            fprintf(stderr, "error finding PASS column\n");
            myexit(status);
        }
        if (fits_read_col(fptr, TINT, colnum, frow, felem, nrows, &nullval,
                          ipass, &anynulls, &status) ) {
            fprintf(stderr, "error reading PASS column\n");
            myexit(status);
        }
        try {
            P.reserve(400000);
        } catch (std::exception & e) {
            myexception(e);
        }
        for (ii = 0; ii < nrows; ii++) {
            //  fprintf(stdout, "in desi %d\n", in_desi[ii]);
            if ( (in_desi[ii] == 1) && (obsconditions[ii] != 0) ) {
                if (ra[ii] <   0.) ra[ii] += 360.;
                if (ra[ii] >= 360.) ra[ii] -= 360.;
                if ( ( dec[ii] < -90.) || ( dec[ii] > 90.) ) {
                    std::cout << "DEC=" << dec << " out of range reading " <<
                    F.tileFile << std::endl;
                    myexit(1);
                }
                double theta = (90.0 - dec[ii]) * M_PI / 180.;
                double phi   = (ra[ii]        ) * M_PI / 180.;
                struct plate Q;
                Q.tileid = tile_id[ii];
                Q.obsconditions = obsconditions[ii];
                // std::cout << "TILEID " << tileid << std::endl;
                Q.tilera        = ra[ii];
                Q.tiledec       = dec[ii];
                Q.nhat[0]    = sin(theta) * cos(phi);
                Q.nhat[1]    = sin(theta) * sin(phi);
                Q.nhat[2]    = cos(theta);
                // <- be careful, format of input file
                Q.ipass      = ipass[ii];
                // <- added
                Q.av_gals.resize(F.Nfiber);
                // <- added
                Q.density.resize(F.Nfiber);
                // was Nfbp
                Q.SS_av_gal.resize(F.Npetal);
                // was Nfbp
                Q.SF_av_gal.resize(F.Npetal);
                Q.SS_in_petal.resize(F.Npetal);
                Q.SF_in_petal.resize(F.Npetal);
                Q.SS_av_gal_fiber.resize(F.Nfiber);
                Q.SF_av_gal_fiber.resize(F.Nfiber);
                for (int i = 0; i < F.Npetal; ++i) {
                    Q.SS_in_petal[i] = 0;
                }
                for (int i = 0; i < F.Npetal; ++i) {
                    Q.SF_in_petal[i] = 0;
                }
                try {
                    P.push_back(Q);
                } catch (std::exception & e) {
                    myexception(e);
                }
            }
        }

        // Free memory
        free(obsconditions);
        free(ipass);
        free(in_desi);
        free(tile_id);
        free(ra);
        free(dec);

        // Close file
        fits_close_file(fptr, &status);

    }
    printf(" size of P  %lu\n", P.size() );
    // Map each valid tileid in order to an index in P[].
    // Tileid is an arbitrary int
    std::map <int, int> invert_tile;
    std::map <int, int>::iterator tileid_to_idx;
    std::pair <std::map <int, int>::iterator, bool> ret;
    init_time_at(time, "# Start invert tiles", t);
    for (unsigned i = 0; i < P.size(); ++i) {
        ret = invert_tile.insert(std::make_pair(P[i].tileid, i) );
        // Check for duplicates (std::map.insert only creates keys, fails on
        // duplicate keys)
        if ( ret.second == false ) {
            std::cerr << "ERROR: Duplicate tileid " << P[i].tileid <<
            " in tileFile!" << std::endl;
            std::exit(1);
            // std::ostringstream o;
            // o << "Duplicate tileid " << P[i].tileid << " in tileFile!";
            // throw std::logic_error(o.str().c_str());
        }
    }
    print_time(time, "# ..inversion  took :");
    // Create PQ, a subset of P containing those tileids specified in the
    // surveyFile, in the order in which they are specified.
    init_time_at(time, "# do inversion of used plates", t);
    for (size_t i = 0; i < survey_list.size(); ++i) {
        tileid        = survey_list[i];
        tileid_to_idx = invert_tile.find(tileid);
        if (tileid_to_idx == invert_tile.end() ) {
            // Can end up with no mapping if surveyFile contains valid tileids
            // that have in_desi = 0 in the tileFile.
            std::cerr << "ERROR: surveyFile contains tileid " << tileid <<
            ", which is not included (or has in_desi = 0) in tileFile." <<
            std::endl;
            std::exit(1);
            // std::ostringstream o;
            // o << "surveyFile contains tileid " << tileid << ", which is not
            // included (or has in_desi = 0) in tileFile.";
            // throw std::range_error(o.str().c_str());
        }
        // Found a valid index, push the tile to the ordered list.
        PQ.push_back(P[tileid_to_idx->second]);
    }
    print_time(time, "# .. sued plates inversion  took :");
    return (PQ);
}

// Assignment
// -----------------------------------------------------------------------------
Assignment::Assignment (const MTL & M, const Feat & F) {
    // galaxy assigned to tile-fiber TF[j][k]
    TF = initTable(F.Nplate, F.Nfiber, -1);
    // tile-fiber pair for galaxy  GL[g]
    GL = initPtable(F.Ngal, 0);
    inv_order = initList(F.Nplate, -1);
    next_plate = 0;
    kinds = initCube(F.Nplate, F.Npetal, F.Categories);
    // initialized to number of fibers on a petal
    unused = initTable(F.Nplate, F.Npetal, F.Nfbp);
}

Assignment::~Assignment () {
}

// Assign g with tile/fiber (j,k), and check for duplicates
void Assignment::assign (int j, int k, int g, MTL & M, Plates & P,
                         const FP & pp) {
    // Assign (j,k)
    int q = TF[j][k];
    if (q != -1) {
        printf(
            "### !!! ### DUPLICATE (j,k) = (%d,%d) assigned with g = %d and %d ---> information on first g lost \n", j, k, q,
            g);
        myexit(1);
    }
    // Assign g
    TF[j][k] = g;
    // pair list, tf's for this g
    Plist pl = GL[g];
    pair p = pair(j, k);
    for (size_t i = 0; i < pl.size(); i++) {
        if (pl[i].f == j) {
            printf(
                "### !!! ### DUPLICATE g = %d assigned with (j,k) = (%d,%d) and (%d,%d) ---> information on first (j,k) lost \n", g,
                pl[i].f, pl[i].s, j, k);
            // Can be commented if want to force continuing
            myexit(1);
        }
    }
    GL[g].push_back(p);
    M[g].nobs_done++;
    M[g].nobs_remain--;
    if (M[g].SF) {
        int q = pp[k].spectrom;
        P[j].SF_in_petal[q] += 1;
    }
    if (M[g].SS) {
        int q = pp[k].spectrom;
        P[j].SS_in_petal[q] += 1;
    }
    unused[j][pp[k].spectrom]--;
}

void Assignment::unassign (int j, int k, int g, MTL & M, Plates & P,
                           const FP & pp) {
    // diagnostic
    if (TF[j][k] == -1) {
        printf(
        "### !!! ### TF (j,k) = (%d,%d) gets unassigned but was already not assigned\n", j, k);
    }
    int a = isfound(pair(j, k), GL[g]);
    if (a == -1) {
        printf(
        "### !!! ### Galaxy g = %d gets unassigned but was already not assigned\n", g);
    }
    TF[j][k] = -1;
    if (a != -1) {
        erase(a, GL[g]);
    }
    M[g].nobs_done--;
    M[g].nobs_remain++;
    if (M[g].SF) {
        int p = pp[k].spectrom;
        P[j].SF_in_petal[p] -= 1;
    }
    if (M[g].SS) {
        int p = pp[k].spectrom;
        P[j].SS_in_petal[p] -= 1;
    }
    unused[j][pp[k].spectrom]++;
}

int Assignment::is_assigned_jg (int j, int g) const {
    // is galaxy g assigned on tile j
    for (size_t i = 0; i < GL[g].size(); i++) {
        if (GL[g][i].f == j) {
            return (int)i;
        }
    }
    return -1;
}

int Assignment::is_assigned_jg (int j, int g, const MTL & M,
                                const Feat & F) const {
    // No occurrence too nearby in tiles
    for (size_t i = 0; i < GL[g].size(); i++) {
        if ( ( fabs(j - GL[g][i].f) < F.InterPlate) ||
             ( j == GL[g][i].f) ) {
            return i;
        }
    }
    return -1;
}

bool Assignment::is_assigned_tf (int j, int k) const {
    return (TF[j][k] != -1);
}

int Assignment::na (const Feat & F, int begin, int size) const {
    // unassigned fibers in tiles begin to begin+size
    int size1 = (size == -1) ? F.Nplate : size;
    int cnt(0);
    for (int j = begin; j < begin + size1; j++) {
        for (int k = 0; k < F.Nfiber; k++) {
            if (TF[j][k] != -1) {
                cnt++;
            }
        }
    }
    return cnt;
}

Plist Assignment::chosen_tfs (int g, const Feat & F, int begin) const {
    // creates list of tile-fibers observing g starting from plate begin
    Plist chosen;
    for (size_t i = 0; i < GL[g].size(); i++) {
        pair tf = GL[g][i];
        if (begin <= tf.f ) {
            if (TF[tf.f][tf.s] != g) {
                printf("ERROR in chosen_tfs\n");
                fl();
            }
            chosen.push_back(tf);
        }
    }
    return chosen;
}

Table Assignment::unused_fbp (const FP & pp, const Feat & F) const {
    // table unused fibers on petal  on tile  j
    Table unused = initTable(F.Nplate, F.Npetal);
    // List Sp = pp.spectrom;
    for (int j = 0; j < F.Nplate; j++) {
        for (int k = 0; k < F.Nfiber; k++) {
            if (!is_assigned_tf(j, k) ) {
                unused[j][pp[k].spectrom]++;
            }
        }
    }
    return unused;
}

// not used
List Assignment::unused_f (const Feat & F) const {
    // total unused fibers
    List unused = initList(F.Nplate);
    for (int j = 0; j < F.Nplate; j++) {
        for (int k = 0; k < F.Nfiber; k++) {
            if (!is_assigned_tf(j, k) ) {
                unused[j]++;
            }
        }
    }
    return unused;
}

int Assignment::unused_f (int j, const Feat & F) const {
    // unused fibers on tile j
    int unused(0);
    for (int k = 0; k < F.Nfiber; k++) {
        if (!is_assigned_tf(j, k) ) {
            unused++;
        }
    }
    return unused;
}

int Assignment::unused_fbp (int j, int k, const FP & pp,
                            const Feat & F) const {
    // unused fibers on petal containing fiber k, tile j
    List fibs = pp.fibers_of_sp[pp[k].spectrom];
    int unused(0);
    for (size_t i = 0; i < fibs.size(); i++) {
        if (!is_assigned_tf(j, fibs[i]) ) {
            unused++;
        }
    }
    return unused;
}

int Assignment::nkind (int j, int k, int kind, const MTL & M, const Plates & P,
                       const FP & pp, const Feat & F, bool pet) const {
    // if pet is false, used petal of fiber k,, if pet is true use petal k
    if (!pet) {
        return kinds[j][pp[k].spectrom][kind];
    } else {
        return kinds[j][k][kind];
    }
}

List Assignment::fibs_unassigned (int j, int pet, const MTL & M, const FP & pp,
                                  const Feat & F) const {
    // list of unassigned fibers on petal pet
    List L;
    List fibs = pp.fibers_of_sp[pet];
    for (int kk = 0; kk < F.Nfbp; kk++) {
        int k = fibs[kk];
        if (!is_assigned_tf(j, k) ) {
            L.push_back(k);
        }
    }
    return L;
}

// Returns the radial distance on the plate (mm) given the angle,
// theta (radians).  This is simply a fit to the data provided.
double plate_dist (const double theta) {
    const double p[4] = {8.297e5, -1750., 1.394e4, 0.0};
    double rr = 0;
    for (int i = 0; i < 4; i++) {
        rr = theta * rr + p[i];
    }
    return rr;
}

// returns the angle (theta) on the plate given the distance
// on the plate (mm)
double plate_angle (double r_plate) {
    double theta;
    double delta_theta = 1E-4;
    double error = 1.0;
    theta = 0.1;
    while (abs(error) > 1E-7) {
        error = plate_dist(theta) - r_plate;
        theta -= (error) /
            ( (plate_dist(theta + delta_theta) - plate_dist(theta) ) /
                delta_theta);
    }
    // fprintf(stdout, "%f %f %f\n", r_plate, theta, plate_dist(theta));
    return theta;
}

// Returns the x-y position on the plate centered at P for galaxy O.
struct onplate change_coords (const struct target & O,
                              const struct plate & P) {
    struct onplate obj;
    // Rotate the "galaxy" vector so that the plate center is at z-hat.
    double nhat1[3], nhat2[3];
    const double ct = P.nhat[2], st = sqrt(1 - P.nhat[2] * P.nhat[2]) + 1e-30;
    const double cp = P.nhat[0] / st, sp = P.nhat[1] / st;
    // First rotate by -(Phi-Pi/2) about z. Note sin(Phi-Pi/2)=-cos(Phi)
    // and cos(Phi-Pi/2)=sin(Phi).
    nhat1[0] =  O.nhat[0] * sp - O.nhat[1] * cp;
    nhat1[1] =  O.nhat[0] * cp + O.nhat[1] * sp;
    nhat1[2] =  O.nhat[2];
    // then rotate by Theta about x
    nhat2[0] =  nhat1[0];
    nhat2[1] =  nhat1[1] * ct - nhat1[2] * st;
    nhat2[2] =  nhat1[1] * st + nhat1[2] * ct;
    // now work out the "radius" on the plate
    double tht = sqrt(sq(nhat2[0], nhat2[1]) );
    double rad = plate_dist(tht);
    // the x-y position is given by our nhat's scaled by this
    obj.pos[0] = nhat2[0] / tht * rad;
    obj.pos[1] = nhat2[1] / tht * rad;
    return obj;
}

struct onplate radec2xy (const struct target & O, const struct plate & P) {
    // following
    // https://github.com/desihub/desimodel/blob/master/py/desimodel/focalplane.py#L259
    struct onplate obj;
    double inc;  // inclination
    double coord[3], coord1[3], coord2[3];
    double deg_to_rad = M_PI / 180.0;
    double x, y, z, x0, y0, z0, x_focalplane, y_focalplane, radius_rad;
    double newteldec, newtelra, ra_rad, dec_rad, q_rad, radius_mm;
    double testra, testdec, dra, ddec;
    double arcsec = 1.0 / 3600.0;
    
    // Inclination is 90 degrees minus the declination in degrees
    inc = 90.0 - O.dec;
    x0 = sin(inc * deg_to_rad) * cos(O.ra * deg_to_rad);
    y0 = sin(inc * deg_to_rad) * sin(O.ra * deg_to_rad);
    z0 = cos(inc * deg_to_rad);
    coord[0] = x0;
    coord[1] = y0;
    coord[2] = z0;
    
    coord1[0] = cos(P.tilera * deg_to_rad) * coord[0] + sin(
        P.tilera * deg_to_rad) * coord[1];
    coord1[1] = -sin(P.tilera * deg_to_rad) * coord[0] + cos(
        P.tilera * deg_to_rad) * coord[1];
    coord1[2] = coord[2];
    
    coord2[0] = cos(P.tiledec * deg_to_rad) * coord1[0] + sin(
        P.tiledec * deg_to_rad) * coord1[2];
    coord2[1] =  coord1[1];
    coord2[2] = -sin(P.tiledec * deg_to_rad) * coord1[0] + cos(
        P.tiledec * deg_to_rad) * coord1[2];
    x = coord2[0];
    y = coord2[1];
    z = coord2[2];
    newteldec = 0;
    newtelra = 0;
    ra_rad = atan2(y, x);
    if (ra_rad < 0) {
        ra_rad = 2.0 * M_PI + ra_rad;
    }
    dec_rad = (M_PI / 2) - acos(z / sqrt( (x * x) + (y * y) + (z * z) ) );
    radius_rad = 2 *
        asin(sqrt( (pow(sin( (dec_rad - newteldec) / 2), 2) ) +
        ( (cos(newteldec) ) * cos(dec_rad) *
        (pow(sin( (   ra_rad - newtelra) / 2), 2) ) ) ) );
    q_rad = atan2(z, -y);
    radius_mm = plate_dist(radius_rad);
    x_focalplane = radius_mm * cos(q_rad);
    y_focalplane = radius_mm * sin(q_rad);
    obj.pos[0] = x_focalplane;
    obj.pos[1] = y_focalplane;
    // test the conversion
    xy2radec(&testra, &testdec, P.tilera, P.tiledec, obj.pos[0], obj.pos[1]);
    dra = (testra * cos(testdec * M_PI / 180.0) - O.ra *
        cos(O.dec * M_PI / 180.0) ) / arcsec;
    ddec = (testdec - O.dec) / arcsec;
    if (fabs(dra) > 0.01) {
        // 0.01 arcsecond precision
        fprintf(stderr,
                "onplate problem with xy2radec conversion [dRA (arcsec)]: %f\n",
                dra);
        fprintf(stderr, "[dDEC (arcsec)]: %f \n", ddec);
        myexit(1);
    }
    if (fabs(ddec) > 0.01) {
        // 0.01 arcsecond precision
        fprintf(stderr,
                "onplate problem with xy2radec conversion [dDEC]: %f\n", ddec);
        fprintf(stderr, "[dRA]: %f\n", dra);
        myexit(1);
    }
    return obj;
}

// Returns the ra-dec position of an x, y position on the focal plane given the
// telescope pointing telra, teldec
void xy2radec (double * ra, double * dec, double telra, double teldec,
               double x, double y) {
    // following
    // https://github.com/desihub/desimodel/blob/master/py/desimodel/focalplane/geometry.py#L165
    double coord[3], coord1[3], coord2[3];
    double x1, y1, z1, x2, y2, z2;
    double teldec_rad = teldec * M_PI / 180.;
    double telra_rad = telra * M_PI / 180.;
    double radius;
    double r_rad;
    double q, q_rad;
    double ra_rad, dec_rad;
    
    // radial distance on the focal plane in radians
    radius = sqrt(x * x + y * y);
    r_rad = plate_angle(radius);
    
    // fprintf(stdout, "tel ra %f tel dec %f\n", telra, teldec);
    // q signifies the angle the position makes with the +x-axis of focal plane
    q_rad = atan2(y, x);
    q = q_rad * 180.0 / M_PI;
    
    // The focal plane is oriented with +yfocal = +dec but +xfocal = -RA
    // Rotate clockwise around z by r_rad
    // zrotate = np.zeros(shape=(3,3))
    // zrotate[0] = [np.cos(r_rad), np.sin(r_rad), 0]
    // zrotate[1] = [-np.sin(r_rad), np.cos(r_rad), 0]
    // zrotate[2] = [0, 0, 1]
    // v1 = zrotate.dot(v0)

    x1 = cos(r_rad);     // y0=0 so drop sin(r_rad) term
    y1 = -sin(r_rad);    // y0=0 so drop cos(r_rad) term
    z1 = 0.0;
    
    // clockwise rotation around the x-axis
    // xrotate = np.zeros(shape=(3,3))
    // q_rad = np.radians(q)
    // xrotate[0] = [1, 0, 0]
    // xrotate[1] = [0, np.cos(q_rad), np.sin(q_rad)]
    // xrotate[2] = [0, -np.sin(q_rad), np.cos(q_rad)]
    
    x2 = x1;
    y2 = y1 * cos(q_rad);           //# z1=0 so drop sin(q_rad) term
    z2 = -y1 * sin(q_rad);    

    coord[0] = x2;
    coord[1] = y2;
    coord[2] = z2;
    
    // Clockwise rotation around y axis by declination of the tile center
    coord1[0] = cos(teldec_rad) * coord[0] - sin(teldec_rad) * coord[2];
    coord1[1] = coord[1];
    coord1[2] = sin(teldec_rad) * coord[0] + cos(teldec_rad) * coord[2];

    // Counter-clockwise rotation around the z-axis by the right ascension of the tile center
    coord2[0] = cos(telra_rad) * coord1[0] - sin(telra_rad) * coord1[1];
    coord2[1] = sin(telra_rad) * coord1[0] + cos(telra_rad) * coord1[1];
    coord2[2] = coord1[2];

    ra_rad = atan2(coord2[1], coord2[0]);
    if (ra_rad < 0) ra_rad = 2.0 * M_PI + ra_rad;
    // fprintf(stdout, "NORM %f %f %f\n", coord4[0], coord4[1], coord4[2]);
    dec_rad = (M_PI / 2.0) - acos(coord2[2]);
    
    *ra = std::fmod( (ra_rad * 180.0 / M_PI), 360.0);
    *dec = dec_rad * 180.0 / M_PI;
    // fprintf(stdout, "FINAL: %f %f \n", *ra, *dec);
}

bool collision (dpair O1, dpair G1, dpair O2, dpair G2, const Feat & F) {
    double dist_sq = sq(G1, G2);
    if (dist_sq < sq(F.Collide) ) return true;
    if (dist_sq > sq(F.NoCollide) ) return false;
    PosP posp(3, 3);
    polygon fh1 = F.fh;
    polygon fh2 = F.fh;
    polygon cb1 = F.cb;
    polygon cb2 = F.cb;
    repos_cb_fh(cb1, fh1, O1, G1, posp);
    repos_cb_fh(cb2, fh2, O2, G2, posp);
    if (collision(fh1, fh2) ) return true;
    if (collision(cb1, fh2) ) return true;
    if (collision(cb2, fh1) ) return true;
    return false;
}

// (On plate p) finds if there is a collision if fiber k would observe galaxy g
// (collision with neighbor)
//  j is in list that runs to F.Nplate since it is used in TF[j][k]
int Assignment::find_collision (int j, int k, int g, const FP & pp,
                                const MTL & M, const Plates & P,
                                const Feat & F, int col) const {
    // check all neighboring fibers
    bool bol = (col == -1) ? F.Collision : false;
    if (bol) {
        return -1;
    }
    dpair G1 = projection(g, j, M, P);
    for (size_t i = 0; i < pp[k].N.size(); i++) {
        // i numbers the fibers neighboring fiber k
        int kn = pp[k].N[i];
        int gn = TF[j][kn];
        if (gn != -1) {
            dpair G2 = projection(gn, j, M, P);
            bool b =
                F.Exact ? collision(pp[k].coords, G1, pp[kn].coords, G2,
                                    F) : (sq(G1, G2) < sq(F.AvCollide) );
            if (b) {
                return kn;
            }
        }
    }
    return -1;
}

bool Assignment::find_collision (int j, int k, int kn, int g, int gn,
                                 const FP & pp, const MTL & M,
                                 const Plates & P, const Feat & F,
                                 int col) const {
    // check two fibers
    bool bol = (col == -1) ? F.Collision : false;
    if (bol) {
        return false;
    }
    dpair G1 = projection(g, j, M, P);
    dpair G2 = projection(gn, j, M, P);
    return F.Exact ? collision(pp[k].coords, G1, pp[k].coords, G2,
                               F) : (sq(G1, G2) < sq(F.AvCollide) );
}

int Assignment::is_collision (int j, int k, const FP & pp, const MTL & M,
                              const Plates & P, const Feat & F) const {
    // find collision for galaxy g
    int g = TF[j][k];
    if (g != -1) {
        return find_collision(j, k, g, pp, M, P, F, 0);
    } else {
        return -1;
    }
}

float Assignment::colrate (const FP & pp, const MTL & M, const Plates & P,
                           const Feat & F, int jend0) const {
    // rate of collisions
    int jend = (jend0 == -1) ? F.Nplate : jend0;
    int col = 0;
    for (int j = 0; j < jend; j++) {
        List done = initList(F.Nfiber);
        for (int k = 0; k < F.Nfiber; k++) {
            if (done[k] == 0) {
                int c = is_collision(j, k, pp, M, P, F);
                if (c != -1) {
                    done[c] = 1;
                    col += 2;
                }
            }
        }
    }
    return percent(col, jend * F.Nfiber);
}

dpair projection (int g, int j, const MTL & M,
                  const Plates & OP) {
    // x and y coordinates for galaxy observed on plate j
    // struct onplate op = change_coords(M[g],OP[j]);
    struct onplate op = radec2xy(M[g], OP[j]);
    return dpair(op.pos[0], op.pos[1]);
}
