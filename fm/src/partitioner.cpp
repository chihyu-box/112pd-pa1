#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include <random>
#include <chrono>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;
using namespace chrono;

void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.emplace_back(new Net(netName));
            _netName2Id.emplace(netName, netId);
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.find(cellName) == _cellName2Id.end()) {
                        int cellId = _cellNum;
                        _cellArray.emplace_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::partition()
{   
    locaInit();
    printSummary();
    gainInit();
    bListInit();
    FM();  
    return;
}

void Partitioner::locaInit() {
    // ===== initialize ===== 
    for(auto cell : _cellArray) {
        if(_maxPinNum < cell->getPinNum()) { // set max pin number
            _maxPinNum = cell->getPinNum();
        }
        if(_cellName2Id[cell->getName()] % 2 == 0) {
            cell->setPart(0); // set cell to A part
            for(auto num : cell->getNetList()) { // increase part count of each net
                _netArray[num]->incPartCount(0);
            }
            _partSize[0]++; // increase size of partition A
        } 
        else {
            cell->setPart(1); // B part
            for(auto num : cell->getNetList()) {
                _netArray[num]->incPartCount(1);
            }
            _partSize[1]++;
        }
    }
    // Set unlocked cell number
    _unlockNum = _cellNum;

    // count cut size
    _cutSize = 0;
    for(auto net : _netArray) {
        if(!(net->getPartCount(0) == 0 || net->getPartCount(1) == 0)) {
            _cutSize++; // count cut size
        }
    }
    return;
}

void Partitioner::gainInit() {
    // Initialize gain
    for(const auto net : _netArray) {
        if(net->getPartCount(0) + net->getPartCount(1) == 1) continue;
        if(net->getPartCount(0) == 0 || net->getPartCount(1) == 0) { // to site = 0
            for(const auto num : net->getCellList()) {
                _cellArray[num]->decGain();
            }
        }
        else {
            if(net->getPartCount(0) == 1) { // from site (A) = 1
                for(const auto num : net->getCellList()) {
                    if(_cellArray[num]->getPart() == 0) {
                        _cellArray[num]->incGain();
                    }
                }
            }
            if(net->getPartCount(1) == 1) { // from site (B) = 1
                for(const auto num : net->getCellList()) {
                    if(_cellArray[num]->getPart() == 1) {
                        _cellArray[num]->incGain();
                    }
                }
            }
        }
    }
    return;
}

void Partitioner::bListInit() {
    // Build bucket list
    for(int i = _maxPinNum; i>=-_maxPinNum; --i) { // Initialize nullptr to bucket list
        Node *headA = new Node(-1);
        headA->setNext(headA); headA->setPrev(headA);
        Node *headB = new Node(-1);
        headB->setNext(headB); headB->setPrev(headB);
        _bList[0].emplace(make_pair(i, headA));
        _bList[1].emplace(make_pair(i, headB));
    }
    for(const auto cell : _cellArray) {
        const auto current = cell->getNode();
        const auto head = _bList[cell->getPart()][cell->getGain()];
        const auto tail = head->getPrev(); // store first node (not count in head)
        insertBList(tail, current, head);
    }
    return;
}

void Partitioner::FM() {
    cout << "start FM" << endl;
    cout << "cutsize " << _cutSize << endl;
    do {
        while(_unlockNum != 0) {
            // cout << "unlocked num " << _unlockNum << endl;
            pickFM();
        } 
        cout << "iter " << ++_iterNum << endl;
        cout << "cutsize " << _cutSize << endl;
        execFM();
        reset();
    } while (_maxAccGain > 0);
    return;
}

inline void Partitioner::pickFM() {
    // pick cell has max gain
    pickMaxGainCell(_maxPinNum);
    auto candiCell = _cellArray[_maxGainCell->getId()];
    // if balanced
    if(!balance(candiCell)) {
        pickPartMaxGainCell(!candiCell->getPart(), _maxPinNum);
        candiCell = _cellArray[_maxGainCell->getId()];
    }
    // record candidate cell for update bucket list
    _moveStack.emplace_back(candiCell->getNode()->getId());
    _gainStack.emplace_back(candiCell->getGain());
    updateCellGain(candiCell);
    moveCell(candiCell);
    lockCell(candiCell);
    updateBList();
    return;
}

inline void Partitioner::execFM() {
    _maxAccGain = 0;
    int bestMoveNum = -1;
    int tempAccGain = 0;
    for(int i=0; i<_gainStack.size(); ++i) {
        tempAccGain += _gainStack[i];
        if(tempAccGain > _maxAccGain) {
            _maxAccGain = tempAccGain;
            bestMoveNum = i;
        }
    }
    if(bestMoveNum >= 0) {
        _moveNum += bestMoveNum+1;
        changedCell.reserve(_moveNum);
        for(int i=0; i<=bestMoveNum; ++i) {
            changedCell.emplace_back(_moveStack[i]);
            updateCellGain(_cellArray[_moveStack[i]]);
            moveCell(_cellArray[_moveStack[i]]);
            _cutSize -= _gainStack[i];
        }
        updateBList();
    }
    return;
}

inline bool Partitioner::balance(Cell* candiCell)
{
    return _partSize[candiCell->getPart()]-1 >= double((1.0-_bFactor)/2.0) * _cellNum &&
           _partSize[!candiCell->getPart()]+1 <= double((1.0+_bFactor)/2.0) * _cellNum;
}

inline void Partitioner::updateCellGain(Cell* candiCell) 
{   
    const auto fromSite = candiCell->getPart();
    const auto toSite = !fromSite;
    const auto candiCellId = candiCell->getNode()->getId(); 
    // emplace candidate cell to changed cell set
    set<int> modifiedCell;
    modifiedCell.emplace(candiCellId);
    // update candidate cell gain
    candiCell->setGain(0-candiCell->getGain());
    // update cells who has connection with candidate cell 
    // Ex. candidate cell is on A part
    for(const auto netId : candiCell->getNetList()) {
        auto fromCount = _netArray[netId]->getPartCount(fromSite);
        auto toCount = _netArray[netId]->getPartCount(toSite);
        if(fromCount + toCount == 0) continue;
        for(const auto cellId : _netArray[netId]->getCellList()) {
            // pass candidate cell, candidate cell gain should be -gain
            if(cellId == candiCellId) continue;
            // store orignial gain
            auto cell = _cellArray[cellId];
            auto cellGain = cell->getGain();
            // before move, to site (B) = 0, from site (A) + 1
            if(toCount == 0) {
                cell->incGain();
            }
            // before move, to site (B) = 1, to site (B) - 1
            else if(toCount == 1 && cell->getPart() == toSite) {
                cell->decGain();
            }
            // after move, from site (A) = 0, to site (B) - 1
            // candidate cell is on from site, and it has not moved to to site yet
            if(fromCount == 1) {
                cell->decGain();
            }
            // after move, from site (A) = 1, from site (A) + 1
            else if(fromCount == 2 && cell->getPart() == fromSite) {
                cell->incGain();
            }
            if(cell->getGain() != cellGain) {
                modifiedCell.emplace(cellId);
            }
        }
    }
    changedCell.insert(changedCell.end(), modifiedCell.begin(), modifiedCell.end());
    return;
}

inline void Partitioner::updateBList() {
    for(const auto cellId : changedCell) {
        // Record prev, next, current node
        auto cell = _cellArray[cellId];
        auto currentNode = cell->getNode();
        auto prevNode = currentNode->getPrev();
        auto nextNode = currentNode->getNext();
        // Connenct previous node to next node  
        prevNode->setNext(nextNode);
        nextNode->setPrev(prevNode);
        // replace cell
        auto headNode = _bList[cell->getPart()][cell->getGain()];
        auto tailNode = headNode->getPrev();
        auto firstNode = headNode->getNext();
        if(cell->getLock()) insertBList(tailNode, currentNode, headNode);
        else insertBList(headNode, currentNode, firstNode);
    }
    changedCell.clear();
    return;
}

inline void Partitioner::moveCell(Cell* candiCell) 
{
    const auto part = candiCell->getPart();
    // update cell info
    candiCell->setPart(!part);
    // update partition info
    _partSize[part]--;
    _partSize[!part]++;
    // update net info
    for(auto netId : candiCell->getNetList()) {
        _netArray[netId]->incPartCount(!part);
        _netArray[netId]->decPartCount(part);
    }
    return;
}

inline void Partitioner::lockCell(Cell *candiCell) {
    candiCell->lock();
    _unlockNum--;
    return;
}

inline void Partitioner::pickMaxGainCell(int gain) {
    for(int i = gain; i>=-_maxPinNum; --i) {
        for (int j = 0; j < 2; ++j) {
            auto headNode = _bList[j][i];
            auto nextNode = headNode->getNext();
            if (nextNode != headNode && !_cellArray[nextNode->getId()]->getLock()) {
                _maxGainCell = nextNode;
                return;
            }
        }
    }
    return;
}

inline void Partitioner::pickPartMaxGainCell(int part, int gain) {
    for(int i = gain; i>=-_maxPinNum; --i) {
        auto headNode = _bList[part][i];
        auto nextNode = headNode->getNext();
        if (nextNode != headNode && !_cellArray[nextNode->getId()]->getLock()) {
            _maxGainCell = nextNode;
            return;
        }
    }
    return;
}

inline void Partitioner::reset() {
    for(auto cell : _cellArray)
        cell->unlock();
    _moveStack.clear();
    _gainStack.clear();
    _unlockNum = _cellNum;
    return;
}

inline void Partitioner::insertBList(Node *head, Node* current, Node* next) {
    head->setNext(current);
    current->setPrev(head);
    current->setNext(next);
    next->setPrev(current);
    return;
}

void Partitioner::printBList() {
    cout << "==================== Bucket List ====================" << endl;
    for(int i=0; i<2; i++) {
        if(i==0) cout << "======================= A part ======================" << endl;
        else     cout << "======================= B part ======================" << endl;
        for(auto b : _bList[i]) {
            cout << setw(2) << b.first << " , head -> ";
            auto node = b.second->getNext();
            while(node != b.second) {
                cout << _cellArray[node->getId()]->getName() << " -> ";
                node = node->getNext();
            }
            if(node->getId() == -1) cout << "head " << endl;
            else cout << "error" << endl;
        }
    }
    cout << "=====================================================" << endl;
    return;
}

void Partitioner::printCellGain() {
    cout << "==================== Cell Gain ====================" << endl;
    for(auto cell : _cellArray) {
        cout << "Cell name : " << cell->getName() << endl;
        if(cell->getLock()) cout << "Cell gain : lock" << endl;
        else cout << "Cell gain : " << cell->getGain() << endl;
    }
    cout << "=====================================================" << endl;
    return;
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
