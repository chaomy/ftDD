/*
 * @Author: chaomy
 * @Date:   2018-04-14 16:08:12
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-11-29 22:25:12
 */

#include "ftDD.hpp"

/*********************************************************
 * write time, stress, density, velocity
 *********************************************************/
void ftDD::writeData(const string& kk) {
  ofstream ofs("chk" + kk + ".txt", std::iostream::out);
  for (int i = 0; i < toltm[kk].size(); i++)
    ofs << std::setprecision(10) << toltm[kk][i] << " " << strss[kk][i] << " "
        << dnstl[kk][i] << " " << vlc[kk][i] << endl;
  ofs.close();
}

/*********************************************************
 * write time, stress, density, velocity and fitted velocity
 *********************************************************/
void ftDD::writeFitData(const vector<string>& ss) {
  writeResults();
  for (const string& kk : ss) {
    ofstream ofs("chk" + kk + ".txt", std::iostream::out);
    for (int i = 0; i < toltm[kk].size(); i++)
      ofs << std::setprecision(10) << toltm[kk][i] << " " << strss[kk][i] << " "
          << dnstl[kk][i] << " " << vlc[kk][i] << " " << ftvlc[kk][i] << endl;
    ofs.close();
  }
}

/*********************************************************
 * write fitting variables
 *********************************************************/
void ftDD::writeResults() {
  std::ofstream ofs("res.txt", std::ofstream::out);
  ofs << "Error = " << dpar.err << " mob = " << dpar.mob
      << " alpha = " << dpar.alp * dpar.alp << " beta = " << dpar.bta * dpar.bta
      << endl;
  ofs.close();

  std::ofstream of2("cnt.txt", std::ofstream::out);
  of2 << "mob " << dpar.mob << endl;
  of2 << "alp " << dpar.alp << endl;
  of2 << "bta " << dpar.bta << endl;
  of2.close();
}