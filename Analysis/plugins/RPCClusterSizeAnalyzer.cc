// Original Author:  Seungjin Yang
//         Created:  Mon, 08 Jul 2024 05:09:27 GMT

// system include files
#include <memory>
#include <utility>

// user include files
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonRPCHitMatch.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "Validation/MuonHits/interface/MuonHitHelper.h"


#include "TTree.h"


class RPCClusterSizeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RPCClusterSizeAnalyzer(const edm::ParameterSet&);
  ~RPCClusterSizeAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;

  std::pair<const reco::MuonRPCHitMatch*, float> findMuonRPCHitMatch(const reco::MuonChamberMatch&);
  std::pair<edm::OwnVector<RPCRecHit>::const_iterator, bool> findRecHit(const reco::MuonRPCHitMatch*, const RPCRecHitCollection::range&);

  // ----------member data ---------------------------
  const edm::EDGetTokenT<edm::View<reco::Muon> > muon_view_token_;
  const edm::EDGetTokenT<RPCRecHitCollection> rec_hit_collection_token_;


  // tree
  TTree* tree_;

  // branches
  const std::vector<std::string> float_branch_name_vec_ = {
    "muon_pt",
    "muon_eta",
    "muon_phi",
    "muon_dXdZ",
    "muon_dYdZ",
    "muon_dXdZErr",
    "muon_dYdZErr",
  };

  const std::vector<std::string> int_branch_name_vec_ = {
    "region",
    "ring",
    "station",
    "sector",
    "layer",
    "subsector",
    "roll",
  };

  std::map<std::string, float> float_branch_map_;
  std::map<std::string, int> int_branch_map_;

};

RPCClusterSizeAnalyzer::RPCClusterSizeAnalyzer(const edm::ParameterSet& parameter_set)
    : muon_view_token_(consumes<edm::View<reco::Muon> >(parameter_set.getParameter<edm::InputTag>("muon"))),
      rec_hit_collection_token_(consumes<RPCRecHitCollection>(parameter_set.getParameter<edm::InputTag>("rpcRecHit"))),
      tree_(nullptr) {
  usesResource(TFileService::kSharedResource);
}

RPCClusterSizeAnalyzer::~RPCClusterSizeAnalyzer() {
}


void RPCClusterSizeAnalyzer::beginJob() {
  edm::Service<TFileService> file_service;
  tree_ = file_service->make<TTree>("tree", "");

  for (const auto key : float_branch_name_vec_) {
    float_branch_map_.emplace(key, 0);
    tree_->Branch(key.c_str(), &float_branch_map_.at(key));
  }

  for (const auto key : int_branch_name_vec_) {
    int_branch_map_.emplace(key, 0);
    tree_->Branch(key.c_str(), &int_branch_map_.at(key));
  }
}

void RPCClusterSizeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muon", edm::InputTag{"muons"});
  desc.add<edm::InputTag>("rpcRecHit", edm::InputTag{"rpcRecHits"});
  descriptions.addDefault(desc);
}

void RPCClusterSizeAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& event_setup) {
  const edm::View<reco::Muon>* muon_view = nullptr;
  if (auto handle = event.getHandle(muon_view_token_)) {
    muon_view = handle.product();
  } else {
    edm::LogError("") << "failed to get View<Muon>";
    return;
  }

  const RPCRecHitCollection* rec_hit_collection = nullptr;
  if (auto handle = event.getHandle(rec_hit_collection_token_)) {
    rec_hit_collection = handle.product();

  } else {
    edm::LogError("") << "failed to get repcRecHits";
    return;
  }

  //---------------------------------------------------------------------------
  for (const auto& muon : *muon_view) {
    for (const auto& muon_chamber_match : muon.matches()) {
      if (not MuonHitHelper::isRPC(muon_chamber_match.id)) {
        continue;
      }

      const auto [muon_hit_match, dx] = findMuonRPCHitMatch(muon_chamber_match);
      if (muon_hit_match == nullptr) {
        edm::LogInfo("") << "no MuonRPCHitMatch";
        continue;
      }

      const RPCDetId rpc_id{muon_chamber_match.id};

      const auto rec_hit_range = rec_hit_collection->get(rpc_id);
      const auto [hit, hit_found] = findRecHit(muon_hit_match, rec_hit_range);
      if (not hit_found) {
        edm::LogError("") << "no RPCRecHit matched with MuonRPCHitMatch";
        continue;
      }

      float_branch_map_.at("muon_pt") = muon.pt();
      float_branch_map_.at("muon_eta") = muon.eta();
      float_branch_map_.at("muon_phi") = muon.phi();

      float_branch_map_.at("muon_dXdZ") = muon_chamber_match.dXdZ;
      float_branch_map_.at("muon_dYdZ") = muon_chamber_match.dYdZ;
      float_branch_map_.at("muon_dXdZErr") = muon_chamber_match.dXdZErr;
      float_branch_map_.at("muon_dYdZErr") = muon_chamber_match.dYdZErr;

      int_branch_map_.at("region") = rpc_id.region();
      int_branch_map_.at("ring") = rpc_id.ring();
      int_branch_map_.at("station") = rpc_id.station();
      int_branch_map_.at("sector") = rpc_id.sector();
      int_branch_map_.at("layer") = rpc_id.layer();
      int_branch_map_.at("subsector") = rpc_id.subsector();
      int_branch_map_.at("roll") = rpc_id.roll();

      tree_->Fill();

    }    // MuonChamberMatch
  } // Muon
}


// find a MuonRPCHitMatch closest to Muon on the chamber
std::pair<const reco::MuonRPCHitMatch*, float>
RPCClusterSizeAnalyzer::findMuonRPCHitMatch(const reco::MuonChamberMatch& muon_chamber_match) {
  const reco::MuonRPCHitMatch* closest_match = nullptr;
  float min_dx = 1e9;

  const float muon_x = muon_chamber_match.x;

  for (const auto& muon_hit_match : muon_chamber_match.rpcMatches) {
    const float dx = std::abs(muon_hit_match.x - muon_x);
    if (dx < min_dx) {
      closest_match = &muon_hit_match;
      min_dx = dx;
    }
  }

  return std::make_pair(closest_match, min_dx);
}


std::pair<edm::OwnVector<RPCRecHit>::const_iterator, bool>
RPCClusterSizeAnalyzer::findRecHit(const reco::MuonRPCHitMatch* muon_hit_match,
                                   const RPCRecHitCollection::range& rec_hit_range) {

  edm::OwnVector<RPCRecHit>::const_iterator matched_hit = rec_hit_range.second;
  for (auto rec_hit = rec_hit_range.first; rec_hit != rec_hit_range.second; ++rec_hit) {
    if ((muon_hit_match->x == rec_hit->localPosition().x()) and (muon_hit_match->bx == rec_hit->BunchX())) {
      matched_hit = rec_hit;
      break;
    }

  }

  return std::make_pair(matched_hit, matched_hit != rec_hit_range.second);
}

DEFINE_FWK_MODULE(RPCClusterSizeAnalyzer);
