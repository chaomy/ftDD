/*
 * @Author: chaomy
 * @Date:   2018-03-16 14:42:04
 * @Last Modified by:   chaomy
 * @Last Modified time: 2019-01-08 12:43:32
 */

#include "ftDD.hpp"

ftDD::ftDD(int argc, char* argv[]) {
  // default setup
  defaultSetup();

  // get the name of input file
  parseArgs(argc, argv);

  // read parameters from the file
  initParm();

  // read DD data
  initData();

  // run optimization
  (this->*run[sparams["ptype"]])(cases);
}

void ftDD::defaultSetup() {
  calobj["ma"] = &ftDD::errLinearMob;
  calobj["a"] = &ftDD::errLinearMobFixM;
  calobj["mab"] = &ftDD::errLinearPrc;
  calobj["b"] = &ftDD::errLinearMobFixMAlp;

  run["ma"] = &ftDD::optErate;
  run["a"] = &ftDD::optFixMob; 
  run["mab"] = &ftDD::optPrec;
  run["b"] = &ftDD::optFixMobAlpha;
  run["0"] = &ftDD::optFixAll;

  sparams["ptype"] = "ma";
  dparams["burg"] = 2.86e-10;
  dparams["taylor"] = 0.40824829046386307;
  dparams["mu"] = 2.7e10;
  dparams["mob"] = 900, dparams["alp"] = 0.1, dparams["bta"] = 10;
}

void ftDD::readData(const string& key) {
  cout << "key is " << key << endl;
  // load flux
  loadDDflx(key);

  // load density
  loadDDdns(key);

  // load property
  loadDDprp(key);

  // times stress by taylor coeff to deresolved stress
  for (auto& ee : strss[key]) ee *= dparams["taylor"];

  for (int i = 0; i < dnstl[key].size(); i++)
    vlc[key].push_back(flxtl[key][i] / (dnstl[key][i] * dparams["burg"]));

  cout << dnstl[key].size() << " " << flxtl[key].size() << " "
       << toltm[key].size() << " " << dnstl[key].size() - flxtl[key].size()
       << endl;

  assert(dnstl[key].size() == flxtl[key].size() &&
         flxtl[key].size() == toltm[key].size());
}

void ftDD::initData() {  // precipitate density
  // prep precipitate data
  prepPrecDat();

  for (const string& key : cases) {
    // read DD data
    readData(key);

    // preprocess DD data
    prepData(key);
  }
}

template <class data_container>
void chop_inputs(data_container& data, int start_point) {
  move(data.begin() + start_point, data.end(), data.begin());
  data.resize(data.size() - start_point);
}

void ftDD::prepData(const string& key) {
  // start fitting from the ns data points (only do flow stress)
  int start_point = 1200;  // ns = 800

  cout << "Total data num " << toltm[key].size() << " use last "
       << toltm[key].size() - start_point << endl;

  chop_inputs(toltm[key], start_point);
  chop_inputs(vlc[key], start_point);
  chop_inputs(strss[key], start_point);
  chop_inputs(dnstl[key], start_point);

  // add weighs since the velocity is different for different cases
  transform(vlc[key].begin(), vlc[key].end(), std::back_inserter(weigh[key]),
            [](double value) { return 1.0 / (1.0 + value); });

  // assign stress to be a constant (average stress)
  double ave_stress =
      std::accumulate(strss[key].begin(), strss[key].end(), 0.0) /
      double(strss[key].size());

  std::fill(strss[key].begin(), strss[key].end(), ave_stress);
}