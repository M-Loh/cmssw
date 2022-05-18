#ifndef DataFormats_HGCalReco_TICLGraph_h
#define DataFormats_HGCalReco_TICLGraph_h

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"

class Node {
    public:
    Node() = default;
    Node(unsigned index, bool isTrackster = true) : index(index), isTrackster(isTrackster) {};
    void addInner(Node &node) {
        innerNodes.push_back(node);
    }
    void addOuter(Node &node) {
        outerNodes.push_back(node);
    }
    ~Node() = default;

    private:
    unsigned index;
    bool isTrackster;
    std::vector<Node> innerNodes;
    std::vector<Node> outerNodes;

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