#ifndef LAGRANGIAN_CHAIN_H
#define LAGRANGIAN_CHAIN_H
#pragma once

#include <lagrangian_object.h>
#include <geometry.h>
#include <cassert>
#include <space.h>

#define REFINE_SORT false // else REFINE_MAX

namespace lagrangian {

class ChainData {
    public:
        template<class V>
        static ChainData createQuad(const V& points, const double& d0, const double& c0, const size_t& dimSpace);
        inline ChainData();
    public: // TODO could be interesting to check c++20 ranges functionalities
        // Links Observers
        inline size_t getLinkNb() const;
        inline size_t getLinkA(const size_t& index) const;
        inline size_t getLinkB(const size_t& index) const;
        // Link setters
        inline void setLinkA(const size_t& index, const size_t& a);
        inline void setLinkB(const size_t& index, const size_t& b);
    public:
        std::vector<std::size_t> links; // contains couples of indexes that points to the corresponding points
        std::vector<std::size_t> cells; // contains group of dimSpace indexs that points to the corresponding links
        std::vector<std::vector<double>> cellArrays;
};

template<class V, template<class> class S> // Vector type, Solver type
class Chain : public Object<V, S> {
    public:
        Chain();
        // Methods
        void initialize();
        void update() override;
        V fState(const V& state, const double& t) const override;
        // Points observers
        size_t getPointNb() const;
        void setX(const size_t& index, const V& x);
        V getX(const size_t& index) const;
        std::vector<size_t> getPointLinks(const size_t& index) const;
        std::vector<size_t> getPointCells(const size_t& index) const;
        // Points setters
        void erasePoint(const size_t& index);
        // Link observers
        double getLinkLength(const size_t& index) const;
        std::vector<size_t> getLinkCells(const size_t& index) const;
        // Link setters
        void eraseLink(const size_t& index);
        // Cell observers
        size_t getCellNb() const;
        std::vector<size_t> getCellLinks(const size_t& index) const;
        std::vector<size_t> getCellPoints(const size_t& index) const;
        V getCellCenter(const size_t& index) const;
        // Cell setters
        void setCellLink(const size_t& index, const size_t& linkIndex, const size_t& link);
        void eraseCell(const size_t& index);
    public:
        // Divide
        virtual void divide(const size_t& c, const size_t& newCell);
        // Compute
        void refine();
        void truncate();
        void repartition();
    public:
        std::shared_ptr<fluid::Fluid<V>> sFluid;
        // Parameters
        // TODO : to remove : double c0;
        double dMin;
        double qMin;
        // Dimension parameters
        size_t dimSpace;
        size_t dimChain;
        size_t nbLinksPerCell;
        // internal data
        ChainData data;
        // Space partition
        space::Bin<size_t, V> spacePartition;
    public:
        // Inherited
        using Object<V, S>::sSolver;
        using Object<V, S>::state;
        using Object<V, S>::t;
};

// Template definition

// ChainData class

ChainData::ChainData() {
    
}

template<class V> // Vector type
ChainData ChainData::createQuad(const V& points, const double& d0, const double& c0, const size_t& dimSpace) {
    ChainData shape;
    // Add links and cells
    for(std::size_t i = 0; i < 3; i++) {
        for(std::size_t e = i + 1; e < 3; e++) {
            shape.links.emplace_back(i);
            shape.links.emplace_back(e);
            shape.cells.emplace_back(shape.getLinkNb() - 1);
        }
    }
    shape.cells.emplace_back(shape.getLinkNb() - 1);
    for(std::size_t i = 1; i < 3; i++) {
        shape.links.emplace_back(i);
        shape.links.emplace_back(3);
        shape.cells.emplace_back(shape.getLinkNb() - 1);
    }
    // Setting qs
    shape.cellArrays.emplace_back();
    shape.cellArrays[0].emplace_back(c0 * d0 * geometry::getTriangleSurface(V(std::begin(points), std::begin(points) + dimSpace), V(std::begin(points)+dimSpace, std::begin(points) + 2*dimSpace), V(std::begin(points)+2*dimSpace, std::begin(points) + 3*dimSpace)));
    shape.cellArrays[0].emplace_back(shape.cellArrays[0][0]);
    // Setting taus
    shape.cellArrays.emplace_back();
    shape.cellArrays[1].emplace_back(0.0);
    shape.cellArrays[1].emplace_back(0.0);
    // Setting c0s
    shape.cellArrays.emplace_back();
    shape.cellArrays[2].emplace_back(c0);
    shape.cellArrays[2].emplace_back(c0);
    // Return shape
    return shape;
}

size_t ChainData::getLinkNb() const {
    return links.size()/2;
}

size_t ChainData::getLinkA(const size_t& index) const {
    return links[2*index];
}

size_t ChainData::getLinkB(const size_t& index) const {
    return links[2*index + 1];
}

void ChainData::setLinkA(const size_t& index, const size_t& a) {
    links[2*index] = a;
}

void ChainData::setLinkB(const size_t& index, const size_t& b) {
    links[2*index + 1] = b;
}

// Chain class

template<class V, template<class T> class S> // Vector type, Solver type, Link type, Cell type
Chain<V, S>::Chain() : Object<V, S>() {
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::initialize() {
    switch(dimChain) {
        case 1 :
            nbLinksPerCell = 1;
            break;
        case 2 :
            nbLinksPerCell = 3;
            break;
        case 3 :
            nbLinksPerCell = 6;
            break;
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
size_t Chain<V, S>::getPointNb() const {
    return state.size()/dimSpace;
}

// TODO should return a const span or range view
template<class V, template<class T> class S> // Vector type, Solver type
V Chain<V, S>::getX(const size_t& index) const {
    V x(dimSpace);
    std::copy(std::begin(state) + index*dimSpace, std::begin(state) + (index + 1)*dimSpace, std::begin(x));
    return x;
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::setX(const size_t& index, const V& x) {
    assert(("Chain::setX, the given point should be od dimension dimSpace", x.size() == dimSpace));
    std::copy(std::begin(x), std::end(x), std::begin(state) + index*dimSpace);
}

template<class V, template<class T> class S> // Vector type, Solver type
std::vector<size_t> Chain<V, S>::getPointLinks(const size_t& pointIndex) const {
    std::vector<size_t> links;
    auto foundPoint = std::begin(data.links);
    while(foundPoint != std::end(data.links)) {
        foundPoint = std::find(foundPoint, std::end(data.links), pointIndex);
        if(foundPoint != std::end(data.links)) {
            links.emplace_back(std::distance(std::begin(data.links), foundPoint)/2);
            foundPoint++;
        }
    }
    return links;
}

template<class V, template<class T> class S> // Vector type, Solver type
std::vector<size_t> Chain<V, S>::getPointCells(const size_t& index) const {
    std::vector<size_t> cells;
    for(const size_t& l : getPointLinks(index)) {
        for(const size_t& c : getLinkCells(l)) {
            if(std::find(std::begin(cells), std::end(cells), c) == std::end(cells))
                cells.emplace_back(c);
        }
    }
    return cells;
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::erasePoint(const size_t& index) {
    state.erase(std::begin(state) + index*dimSpace, std::begin(state) + (index + 1)*dimSpace);
    // Update links
    for(size_t l = 0; l < data.getLinkNb(); l++) {
        const size_t a = data.getLinkA(l);
        const size_t b = data.getLinkB(l);
        assert(("Chain::erasePoint: The erased point shouldn't be linked", a != index and b != index));
        if(a > index)
            data.setLinkA(l, a - 1);
        if(b > index)
            data.setLinkB(l, b - 1);
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
double Chain<V, S>::getLinkLength(const size_t& index) const {
    V vect(dimSpace);
    V ax = getX(data.getLinkA(index));
    V bx = getX(data.getLinkB(index));
    std::transform(std::begin(ax), std::end(ax), std::begin(bx), std::begin(vect), std::minus<double>());
    return geometry::norm<V>(vect);
}

template<class V, template<class T> class S> // Vector type, Solver type
std::vector<size_t> Chain<V, S>::getLinkCells(const size_t& linkIndex) const {
    std::vector<size_t> cells;
    auto foundLink = std::begin(data.cells);
    while(foundLink != std::end(data.cells)) {
        foundLink = std::find(foundLink, std::end(data.cells), linkIndex);
        if(foundLink != std::end(data.cells)) {
            cells.emplace_back(std::distance(std::begin(data.cells), foundLink)/nbLinksPerCell);
            foundLink++;
        }
    }
    return cells;
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::eraseLink(const size_t& index) {
    const size_t& a = data.getLinkA(index);
    size_t b = data.getLinkB(index);
    data.links.erase(std::begin(data.links) + index*2, std::begin(data.links) + (index + 1)*2);
    if(getPointLinks(a).empty()) {
        erasePoint(a);
        if(a < b)
            b--;
    }
    if(getPointLinks(b).empty())
        erasePoint(b);
    // Update cells
    for(size_t c = 0; c < getCellNb(); c++) {
        std::vector<size_t> cLinks = getCellLinks(c);
        for(size_t i = 0; i < cLinks.size(); i++) {
            size_t l = cLinks[i];
            assert(("Chain::eraseLink: The erased link shouldn't be associated to a cell", l != index));
            if(l > index)
                setCellLink(c, i, l - 1);
        }
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
size_t Chain<V, S>::getCellNb() const {
    return data.cells.size()/nbLinksPerCell;
}

// TODO should return a const span or range view
template<class V, template<class T> class S> // Vector type, Solver type
std::vector<size_t> Chain<V, S>::getCellLinks(const size_t& index) const {
    std::vector<size_t> links(nbLinksPerCell);
    std::copy(std::begin(data.cells) + index*nbLinksPerCell, std::begin(data.cells) + (index + 1)*nbLinksPerCell, std::begin(links));
    return links;
}

// TODO should return a const span or range view
template<class V, template<class T> class S>
std::vector<size_t> Chain<V, S>::getCellPoints(const size_t& index) const {
    std::vector<size_t> links = getCellLinks(index);
    std::vector<size_t> points;
    for(size_t l : links) {
        const size_t a = data.getLinkA(l);
        const size_t b = data.getLinkB(l);
        if(std::find(std::begin(points), std::end(points), a) == std::end(points))
            points.emplace_back(a);
        if(std::find(std::begin(points), std::end(points), b) == std::end(points))
            points.emplace_back(b);
    }
    return points;
}

template<class V, template<class T> class S>
V Chain<V, S>::getCellCenter(const size_t& c) const {
    std::vector<size_t> cPoints = getCellPoints(c);
    V cCenter(dimSpace, 0.0);
    cCenter = std::accumulate(std::begin(cPoints), std::end(cPoints), cCenter, [this](V center, const size_t& point) {
        V x = getX(point);
        std::transform(std::begin(x), std::end(x), std::begin(center), std::begin(center), std::plus<double>());
        return center;
    });
    std::transform(std::begin(cCenter), std::end(cCenter), std::begin(cCenter), std::bind(std::multiplies<double>(), std::placeholders::_1, 1/cPoints.size()));
    return cCenter;
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::setCellLink(const size_t& index, const size_t& linkIndex, const size_t& link) {
    data.cells[index*nbLinksPerCell +  linkIndex] = link;
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::eraseCell(const size_t& index) {
    const std::vector<size_t> links = getCellLinks(index);
    data.cells.erase(std::begin(data.cells) + index * nbLinksPerCell, std::begin(data.cells) + (index + 1) * nbLinksPerCell);
    for(auto& array : data.cellArrays) {
        array.erase(std::begin(data.cellArrays[0]) + index);
    }
    
    // TODO, find if
    for(std::size_t l = 0; l < data.getLinkNb(); l++){
        if(getLinkCells(l).empty()) {
            eraseLink(l);
            l--;
        }
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::update() {
    assert(("There should be 1 point or more in the chain", getPointNb() > 0));
    assert(("There should be 1 link or more in the chain", data.getLinkNb() > 0));
    assert(("There should be 1 cell or more in the chain", getCellNb() > 0));
    Object<V, S>::update();
    refine();
    if(qMin > 0.0)
        truncate();
    if(not spacePartition.d.empty())
        repartition();
}

template<class V, template<class T> class S> // Vector type, Solver type
V Chain<V, S>::fState(const V& p_state, const double& p_t) const {
    V dState = p_state;
    for(auto iPoints = std::begin(dState); iPoints != std::end(dState); iPoints += dimSpace){
        const V u = sFluid->getVelocity(V(iPoints, iPoints + dimSpace), p_t);
        std::transform(std::begin(u), std::end(u), iPoints, [](const double& val) { return val; });
    }
    return dState;
}

template<class V, template<class T> class S> // Vector type, Solver type, Link type, Cell type
void Chain<V, S>::refine() {
    // Create links index loop
    std::vector<size_t> links(data.getLinkNb());
    std::iota(std::begin(links), std::end(links), 0);
    // Loop on the links
    size_t biggestLink = *(std::max_element(std::begin(links), std::end(links), [this](const auto& i, const auto& j){
        return getLinkLength(i) < getLinkLength(j);
    }));
    while(getLinkLength(biggestLink) > dMin){
        const size_t a = data.getLinkA(biggestLink);
        const size_t b = data.getLinkB(biggestLink);
        // Add a new point in between TODO interpolation
        const V ax = getX(a);
        const V bx = getX(b);
        std::transform(std::begin(ax), std::end(ax), std::begin(bx), std::back_inserter(state), std::plus<double>());
        std::transform(std::end(state) - dimSpace, std::end(state), std::end(state) - dimSpace, std::bind(std::multiplies<double>(), std::placeholders::_1, 0.5));
        const size_t newPoint = getPointNb() - 1;
        // Add new link
        data.links.emplace_back(newPoint);
        data.links.emplace_back(b);
        const size_t firstNewLink = data.getLinkNb() - 1;
        // Change current link
        data.setLinkB(biggestLink, newPoint);
        // Now updating and adding cells
        const std::vector<size_t> biggestLinkCells = getLinkCells(biggestLink);
        for(const size_t& c : biggestLinkCells) {
            // Add part of new cell
            data.cells.emplace_back(firstNewLink);
            // Update current cell and finish new cell
            const std::vector<size_t> cLinks = getCellLinks(c);
            for(size_t i = 0; i < cLinks.size(); i++){
                const size_t l = cLinks[i];
                if(data.getLinkA(l) == b){
                    // Add new link
                    data.links.emplace_back(newPoint);
                    data.links.emplace_back(data.getLinkB(l));
                    const size_t newLink = data.getLinkNb() - 1;
                    // Add part of new cell
                    data.cells.emplace_back(newLink);
                    data.cells.emplace_back(l);
                    // Update current cell
                    setCellLink(c, i, newLink);
                }
                else if(data.getLinkB(l) == b){
                    // Add new link
                    data.links.emplace_back(data.getLinkA(l));
                    data.links.emplace_back(newPoint);
                    const size_t newLink = data.getLinkNb() - 1;
                    // Add part of new cell
                    data.cells.emplace_back(newLink);
                    data.cells.emplace_back(l);
                    // Update current cell
                    setCellLink(c, i, newLink);
                }
                else if(data.getLinkA(l) != a && data.getLinkB(l) != a){
                    data.cells.emplace_back(l);
                }
            }
            // New cell
            const size_t newCell = getCellNb() - 1;
            // Divide
            divide(c, newCell);
        }
        links = std::vector<size_t>(data.getLinkNb());
        std::iota(std::begin(links), std::end(links), 0);
        biggestLink = *(std::max_element(std::begin(links), std::end(links), [this](const auto& i, const auto& j){
            return getLinkLength(i) < getLinkLength(j);
        }));
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::divide(const size_t& c, const size_t& newCell) {
    // Dividing TODO Should be dependant of the dimChain
    const double q = data.cellArrays[0][c];
    // Getting cells area
    const std::vector<size_t> cPoints = getCellPoints(c);
    const double cArea = geometry::getTriangleSurface(getX(cPoints[0]), getX(cPoints[1]), getX(cPoints[2]));
    const std::vector<size_t> newCellPoints = getCellPoints(newCell);
    const double newCellArea = geometry::getTriangleSurface(getX(newCellPoints[0]), getX(newCellPoints[1]), getX(newCellPoints[2]));
    const double totArea = newCellArea + cArea;
    // Setting new qs
    data.cellArrays[0][c] = q * cArea/totArea;
    data.cellArrays[0].emplace_back(q * newCellArea/totArea);
    for(size_t i = 1; i < data.cellArrays.size(); i++) {
        data.cellArrays[i].emplace_back(data.cellArrays[i][c]);
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::truncate() {
    // Loop on the cells
    for(std::size_t c = 0; c < getCellNb(); c++){
        if(data.cellArrays[0][c] < qMin) { // If quantity is to low, then just remove the cell
            // Removing cell
            eraseCell(c);
            c--;
        }
    }
}

template<class V, template<class T> class S> // Vector type, Solver type
void Chain<V, S>::repartition() {
    spacePartition.clear();
    // Loop on the cells
    for(std::size_t c = 0; c < getCellNb(); c++) {
        spacePartition.add(c, getCellCenter(c));
    }
}

}

#endif
