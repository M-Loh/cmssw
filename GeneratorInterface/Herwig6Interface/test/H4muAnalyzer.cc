/* This is en example for an Analyzer of a Herwig HepMCProduct
   and looks for muons pairs and fills a histogram
   with the invaraint mass of the four. 
*/

//
// Original Author:  Fabian Stoeckli
//         Created:  Tue Nov 14 13:43:02 CET 2006
//
//

// system include files
#include <iostream>
#include <cmath>

#include <TH1D.h>

#include <HepMC/WeightContainer.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//

class H4muAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit H4muAnalyzer(const edm::ParameterSet&);
  ~H4muAnalyzer() override = default;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------

  const edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;
  TH1D* weight_histo;
  TH1D* invmass_histo;
};

H4muAnalyzer::H4muAnalyzer(const edm::ParameterSet& iConfig)
    : hepMCToken_(consumes<edm::HepMCProduct>(edm::InputTag("VtxSmeared"))) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  invmass_histo = fs->make<TH1D>("invmass_histo", "invmass_histo", 60, 170, 180);
}

// ------------ method called to for each event  ------------
void H4muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get HepMC::GenEvent ...
  const edm::Handle<edm::HepMCProduct>& evt_h = iEvent.getHandle(hepMCToken_);
  const HepMC::GenEvent* evt = evt_h->GetEvent();

  // look for stable muons
  std::vector<HepMC::GenParticle*> muons;
  for (auto it = evt->particles_begin(); it != evt->particles_end(); ++it) {
    if (std::abs((*it)->pdg_id()) == 13 && (*it)->status() == 1)
      muons.push_back(*it);
  }

  // if there are at least four muons
  // calculate invarant mass of first two and fill it into histogram
  math::XYZTLorentzVector tot_momentum;
  double inv_mass = 0.0;
  if (muons.size() > 3) {
    for (unsigned int i = 0; i < 4; ++i)
      tot_momentum += muons[i]->momentum();
    inv_mass = tot_momentum.mass();
  }
  invmass_histo->Fill(inv_mass);
}

//define this as a plug-in
typedef H4muAnalyzer H4muExampleAnalyzer;
DEFINE_FWK_MODULE(H4muExampleAnalyzer);
