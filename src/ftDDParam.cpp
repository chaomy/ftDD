/*
 * @Author: chaomy
 * @Date:   2018-06-21 21:09:38
 * @Last Modified by:   chaomy
 * @Last Modified time: 2019-01-08 11:44:01
 */

#include "ftDD.hpp"

/*********************************************************
 * read two type of files :
 * 1. materials paramters
 * 2. infile that specifies job types
 *********************************************************/
void ftDD::readParam() {
  cout << "parfile " << sparams["parfile"] << endl;
  cout << "infile " << sparams["infile"] << endl;

  // read parameter file : parameter file starts with ft.xxx
  ifstream fid(sparams["parfile"], std::ifstream::in);
  vector<string> segs;
  string buff;

  while (getline(fid, buff)) {
    segs.clear();
    split(buff, " ", segs);
    if (sparams.find(segs[0]) != sparams.end()) {
      sparams[segs[0]] = segs[1];
    } else if (dparams.find(segs[0]) != dparams.end()) {
      dparams[segs[0]] = stof(segs[1]);
    }
  }
  fid.close();

  // read cases : the case file start with in.xxx
  ifstream ifs(sparams["infile"], std::ifstream::in);
  getline(ifs, buff);
  std::stringstream s(move(buff));
  while (s >> buff) cases.push_back(buff);
  ifs.close();
}

/*********************************************************
 * Initilize materials parameters
 *********************************************************/
void ftDD::initParm() {
  readParam();

  // copy materials constants to defined object since hashing is slow
  dpar.mob = dparams["mob"], dpar.alp = dparams["alp"],
  dpar.bta = dparams["bta"];

  mcns.taylor = dparams["taylor"], mcns.mu = dparams["mu"];
  mcns.burg = dparams["burg"], mcns.mub = dparams["mu"] * dparams["burg"];

  cout << "Init mob " << dpar.mob << " alp " << dpar.alp << " bta " << dpar.bta
       << endl;
}

/*********************************************************
 * parse arguments to take file names
 *********************************************************/
void ftDD::parseArgs(int argc, char* argv[]) {
  for (int i = 0; i < argc; i++) {
    // parameter file : specify materials constants
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    // input file : specify job to run
    if (!strcmp(argv[i], "--i") || !strcmp(argv[i], "-i"))
      sparams["infile"] = string(argv[++i]);
  }
}
