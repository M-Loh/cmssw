#ifndef DataFormats_HGCalReco_TICLGraph_h
#define DataFormats_HGCalReco_TICLGraph_h

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"

class Node {
    public:
    Node() = default;
    Node(unsigned index, bool isTrackster) : index(index), isTrackster(isTrackster) {};
    void addLinked(Node &node) {
        linkedNodes.push_back(node);
    }
    ~Node() = default;

    private:
    unsigned index;
    bool isTrackster;
    std::vector<Node> linkedNodes;

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