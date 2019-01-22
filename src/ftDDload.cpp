/*
 * @Author: chaomingyang
 * @Date:   2018-03-15 04:06:08
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-12-28 19:06:02
 */

#include "ftDD.hpp"
#define SMALL 1e-12

/*********************************************************
 * load ParaDiS flux results (FCC)
 *********************************************************/
void ftDD::loadDDflx(const string& kk) {  // load dislocation flux
  vector<string> fnms({"fluxtot_b1", "fluxtot_b2", "fluxtot_b3", "fluxtot_b4",
                       "fluxtot_b5", "fluxtot_b6"});
  string buff;
  // each file stores info of one burger, read each of them
  for (int j = 0; j < fnms.size(); j++) {
    ifstream ifs((kk + "_results/fluxdata/" + fnms[j]).c_str(),
                 std::iostream::in);
    double tm[8];
    while (getline(ifs, buff)) {
      sscanf(buff.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf", &tm[0], &tm[1],
             &tm[2], &tm[3], &tm[4], &tm[5], &tm[6], &tm[7]);
      for (int i : {0, 1, 2, 3})
        dd[kk].slips[4 * j + i].flx.push_back(fabs(tm[3 + i]));
    }
    ifs.close();
  }

  // burger vector outer products normal vector
  vector<vector<double>> BON(3, vector<double>(3, 0));

  // normal vector outer products burger vector
  vector<vector<double>> NOB(3, vector<double>(3, 0));

  // sum up those dislocation that have been activated by
  // calculate the resolved matrix 0.5(bxn + nxb)
  flxtl[kk] = vector<double>(dd[kk].slips[0].flx.size());
  for (slip& ee : dd[kk].slips) {
    outProd(ee.b, ee.n, BON);
    outProd(ee.n, ee.b, NOB);
    if (abs(BON[2][2] + NOB[2][2]) > SMALL)
      for (int i = 0; i < flxtl[kk].size(); i++) flxtl[kk][i] += ee.flx[i];
  }
}

/*********************************************************
 * load ParaDiS dislocation density results (FCC)
1         Plastic strain
2         Strain
3         Dislocation density
4         Deleted dislocation density
5         Average dislocation velocity
6         Std deviation of dislocation velocity
7         File version number
8         Burgers vectors [1 1 0] [-1 -1 0]
9         Burgers vectors [-1 1 0] [1 -1 0]
10        Burgers vectors [1 0 1] [-1 0 -1]
11        Burgers vectors [-1 0 1] [1 0 -1]
12        Burgers vectors [0 1 1] [0 -1 -1]
13        Burgers vectors [0 -1 1] [0 1 -1]
14        All other burgers vectors
 *********************************************************/
void ftDD::loadDDdns(const string& kk) {
  string buff;
  ifstream ifs(kk + "_results/properties/density", std::iostream::in);
  double tm[14];
  while (getline(ifs, buff)) {
    sscanf(buff.c_str(),
           "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tm[0],
           &tm[1], &tm[2], &tm[3], &tm[4], &tm[5], &tm[6], &tm[7], &tm[8],
           &tm[9], &tm[10], &tm[11], &tm[12], &tm[13]);
    for (int i : {0, 1, 2, 3, 4, 5})
      dd[kk].slips[4 * i].dns.push_back(fabs(tm[7 + i]));
  }
  ifs.close();

  // burger vector outer products normal vector
  vector<vector<double>> BON(3, vector<double>(3, 0));

  // normal vector outer products burger vector
  vector<vector<double>> NOB(3, vector<double>(3, 0));

  dnstl[kk] = vector<double>(dd[kk].slips[0].dns.size());
  for (int i : {0, 1, 2, 3, 4, 5}) {
    slip& ee = dd[kk].slips[4 * i];
    outProd(ee.b, ee.n, BON);
    outProd(ee.n, ee.b, NOB);
    if (abs(BON[2][2] + NOB[2][2]) > SMALL)
      for (int i = 0; i < dnstl[kk].size(); i++) dnstl[kk][i] += ee.dns[i];
  }
}

/*********************************************************
 * load ParaDiS stress-strain results
 # 1         Strain
 # 2         Stress
 # 3         Elapsed simulation time
 # 4         Curent simulation cycle number
 *********************************************************/
void ftDD::loadDDprp(const string& kk) {  // load strain
  string buff;
  ifstream ifs(kk + "_results/properties/stress_Plastic_strain",
               std::iostream::in);
  double tm[3];
  while (getline(ifs, buff)) {
    sscanf(buff.c_str(), "%lf %lf", &tm[0], &tm[1]);
    strss[kk].push_back(tm[1]);
  }
  ifs.close();
  ifs.open(kk + "_results/properties/time_Plastic_strain", std::iostream::in);
  while (getline(ifs, buff)) {
    sscanf(buff.c_str(), "%lf %lf", &tm[0], &tm[1]);
    toltm[kk].push_back(tm[0]);
  }
  ifs.close();
}
