import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring( 'file:step3.root' )
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('TICLGraph.root')
                                   )

process.demo = cms.EDAnalyzer('DemoAnalyzer',
   tracks    = cms.untracked.InputTag('generalTracks'),
                              )

process.p = cms.Path(process.demo)
