
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace ticl;
using reco::TrackCollection;

class TICLGraphAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  typedef math::XYZVector Vector;
  explicit TICLGraphAnalyzer(const edm::ParameterSet&);
  ~TICLGraphAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  const edm::InputTag trackstersClue3d_;
  const edm::InputTag tracksterGraph_;

  edm::EDGetTokenT<std::vector<Trackster>> tracksters_clue3d_token_;
  edm::EDGetTokenT<TICLGraph> trackster_graph_token_;

  TTree *tree;

  std::vector<math::XYZVector> barycenter;
  std::vector<float> trackster_energy;
  std::vector<std::vector<unsigned int>> node_adj;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TICLGraphAnalyzer::TICLGraphAnalyzer(const edm::ParameterSet& iConfig)
    : trackstersClue3d_(iConfig.getParameter<edm::InputTag>("trackstersclue3d")),
      tracksterGraph_(iConfig.getParameter<edm::InputTag>("tracksterGraph")) {
  edm::ConsumesCollector&& iC = consumesCollector();
  tracksters_clue3d_token_ = iC.consumes<std::vector<ticl::Trackster>>(trackstersClue3d_);
  trackster_graph_token_ = iC.consumes<TICLGraph>(tracksterGraph_);

  edm::Service<TFileService> fs;
  
  auto tree = fs->make<TTree>("tree","Title");

  tree->Branch("barycenter", &barycenter);
  tree->Branch("raw_energy", &trackster_energy);
  tree->Branch("adj", &node_adj);
}

TICLGraphAnalyzer::~TICLGraphAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void TICLGraphAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<std::vector<Trackster>> trackstersclue3d_h;
  iEvent.getByToken(tracksters_clue3d_token_, trackstersclue3d_h);
  auto trackstersclue3d = *trackstersclue3d_h;

  edm::Handle<TICLGraph> trackster_graph_h;
  iEvent.getByToken(trackster_graph_token_, trackster_graph_h);
  auto trackster_graph = *trackster_graph_h;

  auto graph_nodes = trackster_graph.getNodes();

  barycenter.clear();
  trackster_energy.clear();
  node_adj.clear();

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];
    barycenter.push_back(t.barycenter());
    trackster_energy.push_back(t.raw_energy());

    auto node = trackster_graph.getNode(id_t);
    node_adj.push_back(node.getOuter());
  }


  // fill the adjecency list


  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void TICLGraphAnalyzer::beginJob() {
  // please remove this method if not needed


}

// ------------ method called once each job just after ending the event loop  ------------
void TICLGraphAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TICLGraphAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TICLGraphAnalyzer);
