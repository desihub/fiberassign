#ifndef FIBERASSIGN_H
#define FIBERASSIGN_H

#include <string>


void fiberassign_exec (
    std::string mtlfile,
    std::string skyfile,
    std::string stdstarfile,
    std::string surveyfile,
    std::string tilefile,
    std::string fiberfile,
    std::string statusfile,
    std::string rundate,
    std::string outdir,
    long starmask
);


#endif
