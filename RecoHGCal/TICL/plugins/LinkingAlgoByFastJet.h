

#include <memory>
#include <array>
#include "RecoHGCal/TICL/plugins/LinkingAlgoBase.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Math/interface/Vector3D.h"
//#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
//
// class declaration
//
namespace ticl{
  class LinkingAlgoByFastJet final : public LinkingAlgoBase {
  public:
    LinkingAlgoByFastJet(const edm::ParameterSet&);
    ~LinkingAlgoByFastJet() override;

    void linkTracksters(const edm::Handle<std::vector<reco::Track>>,
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const std::vector<reco::Muon> &,
                        const edm::Handle<std::vector<Trackster>>,
                        std::vector<TICLCandidate> &);

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

  private:
    typedef math::XYZVector Vector;

    const double antikt_radius_;
    const double pid_threshold_;
    const double energy_em_over_total_threshold_;
    const std::vector<int> filter_on_categories_;
  };
}

