// -*- C++ -*-
//
// Package:    link/LinkingAlgoByFastJet
// Class:      LinkingAlgoByFastJet
//
/**\class LinkingAlgoByFastJet LinkingAlgoByFastJet.cc RecoHGCal/TICL/plugins/LinkingAlgoByFastJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maximilian Leon Lohmann
//         Created:  Mon, 08 Aug 2022 07:30:37 GMT
//
//

// system include files
#include <memory>

#include <memory>
#include <array>
#include "RecoHGCal/TICL/plugins/LinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/LinkingAlgoByFastJet.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace ticl;

//
// constructors and destructor
//
LinkingAlgoByFastJet::LinkingAlgoByFastJet(const edm::ParameterSet& conf)
  : LinkingAlgoBase(conf),
    filter_on_categories_(conf.getParameter<std::vector<int>>("filter_hadronic_on_categories")),
    pid_threshold_(conf.getParameter<double>("pid_threshold")),
    energy_em_over_total_threshold_(conf.getParameter<double>("energy_em_over_total_threshold")),
    antikt_radius_(conf.getParameter<double>("antikt_radius")) {};


LinkingAlgoByFastJet::~LinkingAlgoByFastJet() {

}

//
// member functions 
//

void LinkingAlgoByFastJet::linkTracksters(const edm::Handle<std::vector<reco::Track>> tkH,
                                          const edm::ValueMap<float> &tkTime,
                                          const edm::ValueMap<float> &tkTimeErr,
                                          const edm::ValueMap<float> &tkTimeQual,
                                          const std::vector<reco::Muon> &muons,
                                          const edm::Handle<std::vector<Trackster>> tsH,
                                          std::vector<TICLCandidate> &resultLinked) {

  constexpr double mpion = 0.13957;
  constexpr float mpion2 = mpion * mpion;

  // const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto isHadron = [&](const Trackster &t) -> bool {
    auto cumulative_prob = 0.;
    for (const auto index : filter_on_categories_) {
      cumulative_prob += t.id_probabilities(index);
    }
    return ((cumulative_prob <= pid_threshold_) and (t.raw_em_energy() != t.raw_energy())) or
           (t.raw_em_energy() < energy_em_over_total_threshold_ * t.raw_energy());
  };

  std::vector<fastjet::PseudoJet> fjInputs_ts;

  for (auto ts : tracksters){
    auto direction = ts.barycenter().Unit();
    direction *= ts.raw_energy();
    auto fpj = fastjet::PseudoJet(direction.X(), direction.Y(), direction.Z(), ts.raw_energy());
    fjInputs_ts.push_back(fpj);
  }

  fastjet::ClusterSequence sequence(fjInputs_ts, JetDefinition(antikt_algorithm, antikt_radius_));
  auto jets_ts = fastjet::sorted_by_pt(sequence.inclusive_jets(0));



  // std::vector<fastjet::PseudoJet> fjInputs_tk;

  // for (auto tk : tracks){
  //   auto direction = tk.barycenter().Unit();
  //   direction *= tk.raw_energy();
  //   auto fpj = fastjet::PseudoJet(direction.X(), direction.Y(), direction.Z(), tk.raw_energy());
  //   fjInputs_tk.push_back(fpj);
  // }
  
  // fastjet::ClusterSequence sequence(fjInputs_tk, JetDefinition(antikt_algorithm, antikt_radius_));
  // auto jets_tk = fastjet::sorted_by_pt(sequence.inclusive_jets(0));


  //create neutral candidates for all jets of tracksters
  std::vector<TICLCandidate> neutralCandidates;
  //neutralMask
  for (const auto &jts : jets_ts) {

    TICLCandidate neutralCandidate;

    //single trackster in a jet
    if (jts.constituents().size() == 1){
      TICLCandidate neutralNoLinks;
      auto ts_i = jts.constituents()[0].user_index();
      neutralNoLinks.addTrackster(edm::Ptr<Trackster>(tsH, ts_i));
      neutralCandidates.push_back(neutralNoLinks);
    }

    else if (jts.constituents().size() > 1){
      for (const auto &component : jts.constituents()){
        neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, component.user_index()));
      }
    }
  }

  for (auto &cand : neutralCandidates) {
    bool isHAD = false;
    double rawE = 0.;
    const auto track = cand.trackPtr();
    double wtSum_baryc[3] = {0};
    for (const auto ts : cand.tracksters()) {
      if (isHadron(*ts))
        isHAD = true;
      rawE += ts->raw_energy();
      wtSum_baryc[0] += (ts->raw_energy()) * (ts->barycenter().x());
      wtSum_baryc[1] += (ts->raw_energy()) * (ts->barycenter().y());
      wtSum_baryc[2] += (ts->raw_energy()) * (ts->barycenter().z());
    }
    Vector combined_baryc(wtSum_baryc[0] / rawE, wtSum_baryc[1] / rawE, wtSum_baryc[2] / rawE);

    if (isHAD) {  // neutral hadron
      cand.setCharge(0);
      cand.setPdgId(130);
      cand.setRawEnergy(rawE);
      float momentum = std::sqrt(rawE * rawE - mpion2);
      math::XYZTLorentzVector p4(momentum * combined_baryc.unit().x(),
                                 momentum * combined_baryc.unit().y(),
                                 momentum * combined_baryc.unit().z(),
                                 rawE);
      cand.setP4(p4);
    } else {  // photon
      cand.setCharge(0);
      cand.setPdgId(22);
      cand.setRawEnergy(rawE);
      math::XYZTLorentzVector p4(
          rawE * combined_baryc.unit().x(), rawE * combined_baryc.unit().y(), rawE * combined_baryc.unit().z(), rawE);
      cand.setP4(p4);
    }
  }

  resultLinked.insert(std::end(resultLinked), std::begin(neutralCandidates), std::end(neutralCandidates));

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void LinkingAlgoByFastJet::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<double>("pid_threshold", 0.5);
  desc.add<double>("energy_em_over_total_threshold", 0.9);
  desc.add<double>("antikt_radius", 0.09)->setComment("Radius to be used while running the Anti-kt clustering");
  desc.add<std::vector<int>>("filter_hadronic_on_categories", {0, 1});

  LinkingAlgoBase::fillPSetDescription(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LinkingAlgoByFastJet);
