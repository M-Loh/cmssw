#ifndef DataFormats_HGCalReco_TICLGraph_h
#define DataFormats_HGCalReco_TICLGraph_h

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"

class Node {
    public:
    Node() = default;
    Node(unsigned index, bool isTrackster = true) : index(index), isTrackster(isTrackster) {};
    void addInner(unsigned int trackster_id) {
        innerNodes.push_back(trackster_id);
    }
    void addOuter(unsigned int trackster_id) {
        outerNodes.push_back(trackster_id);
    }
    ~Node() = default;

    private:
    unsigned index;
    bool isTrackster;
    std::vector<unsigned int> innerNodes;
    std::vector<unsigned int> outerNodes;

};

class TICLGraph {
    public:
    TICLGraph() = default;
    TICLGraph(std::vector<Node> &n) {
        nodes = n;
    };
    ~TICLGraph() = default;

    private:
    std::vector<Node> nodes;

};

#endif