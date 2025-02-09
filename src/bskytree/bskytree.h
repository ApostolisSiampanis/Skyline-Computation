/*
 * bskytree.h
 *
 *  Created on: Jul 6, 2014
 *      Author: dariuss
 *
 */

#ifndef BSKYTREE_H_
#define BSKYTREE_H_

#include <map>
#include <vector>

#include <common/skyline_i.h>
#include <bskytree/node.h>

using namespace std;

class BSkyTree final : public SkylineI {
public:
    BSkyTree(uint32_t num_threads, uint32_t n, uint32_t d);

    ~BSkyTree() override;

    void Init(float **dataset) override;

    vector<int> Execute() override;

private:
    void BSkyTreeS_ALGO();

    void DoPartitioning();

    // PivotSelection methods
    void SelectBalanced();

    static vector<float> SetRangeList(const vector<float> &min_list,
                                      const vector<float> &max_list);

    static float ComputeDistance(const float *value, const vector<float> &min_list,
                                 const vector<float> &range_list);

    bool EvaluatePoint(uint32_t pos) const;

    static void *parallel_execution(void *arg);

    const uint32_t num_threads_;
    const uint32_t n_;
    const uint32_t d_;
    vector<TUPLE_S> data_;

    vector<int> skyline_;
    vector<int> eqm_; // "equivalence matrix"
};

#endif /* BSKYTREE_H_ */
