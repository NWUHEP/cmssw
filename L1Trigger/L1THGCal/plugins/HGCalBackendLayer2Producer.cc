#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerSums.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"

#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include <memory>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // Bueghly
#include <cmath> // Bueghly

class HGCalBackendLayer2Producer : public edm::stream::EDProducer<> {  
 public:    
  HGCalBackendLayer2Producer(const edm::ParameterSet&);
  ~HGCalBackendLayer2Producer() override { }
  
  void beginRun(const edm::Run&, 
                        const edm::EventSetup&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  // inputs
  edm::EDGetToken input_clusters_;
  edm::EDGetToken gen_particles_; // Bueghly
  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;

  std::unique_ptr<HGCalBackendLayer2ProcessorBase> backendProcess_;
};

DEFINE_FWK_MODULE(HGCalBackendLayer2Producer);

HGCalBackendLayer2Producer::
HGCalBackendLayer2Producer(const edm::ParameterSet& conf): 
  input_clusters_(consumes<l1t::HGCalClusterBxCollection>(conf.getParameter<edm::InputTag>("InputCluster")))
{
  gen_particles_ = consumes<reco::GenParticleCollection>(conf.getParameter<edm::InputTag>("GenParticles")); // Bueghly
  //setup Backend parameters
  const edm::ParameterSet& beParamConfig = conf.getParameterSet("ProcessorParameters");
  const std::string& beProcessorName = beParamConfig.getParameter<std::string>("ProcessorName");
  HGCalBackendLayer2ProcessorBase* beProc = HGCalBackendLayer2Factory::get()->create(beProcessorName, beParamConfig);
  backendProcess_.reset(beProc);

  produces<l1t::HGCalMulticlusterBxCollection>(backendProcess_->name());
}

void HGCalBackendLayer2Producer::beginRun(const edm::Run& /*run*/,
                                          const edm::EventSetup& es) 
{                 
  es.get<CaloGeometryRecord>().get("",triggerGeometry_);
  backendProcess_->setGeometry(triggerGeometry_.product());
}

void HGCalBackendLayer2Producer::produce(edm::Event& e, const edm::EventSetup& es) 
{
  // Output collections
  auto be_multicluster_output = std::make_unique<l1t::HGCalMulticlusterBxCollection>();

  // Input collections   
  edm::Handle<l1t::HGCalClusterBxCollection> trigCluster2DBxColl;

  e.getByToken(input_clusters_, trigCluster2DBxColl);

  backendProcess_->run(trigCluster2DBxColl, *be_multicluster_output, es);

  // Block by Bueghly
  bool gen_match_mcls = true;
  if (gen_match_mcls) {
      auto be_multicluster_output_slimmed = std::make_unique<l1t::HGCalMulticlusterBxCollection>();
      edm::Handle<reco::GenParticleCollection> genParticleColl; 
      e.getByToken(gen_particles_, genParticleColl);
      float max_gen_mcl_dr = 0.6;
      for (const auto& particle : *genParticleColl) {
          //cout << "particle pdgid " << particle.pdgId() << endl;
          //cout << "particle status " << particle.status() << endl;
          float gen_eta = particle.eta();
          float gen_phi = particle.phi();
          for (const auto& multicluster : *be_multicluster_output) {
              float mcl_eta = multicluster.eta();
              float mcl_phi = multicluster.phi();
              float deta = fabs(gen_eta - mcl_eta);
              float dphi = fabs(gen_phi - mcl_phi);
              if (dphi > M_PI) 
                  dphi = 2*M_PI - dphi;
              float dr = sqrt(pow(deta, 2) + pow(dphi, 2));
              if (dr < max_gen_mcl_dr) 
                  be_multicluster_output_slimmed.get()->push_back(0, multicluster);
          }             
      }
      e.put(std::move(be_multicluster_output_slimmed), backendProcess_->name()); // Bueghly
  }
  else 
      e.put(std::move(be_multicluster_output), backendProcess_->name());  
}
