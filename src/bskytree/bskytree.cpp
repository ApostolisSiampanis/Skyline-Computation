/*
 * bskytree.cpp
 *
 *  Created on: Jul 6, 2014
 *      Author: dariuss
 */

#include <bskytree/bskytree.h>

#include <algorithm>
#include <cassert>

#include "common/common.h"
#include "util/timing.h"

pthread_barrier_t barrier;
uint32_t head, tail;


struct ThreadData {
    vector<TUPLE_S> *S;
    volatile bool *dead;
    vector<int> *skyline;
    uint32_t thread_id, num_of_threads;

    ThreadData() = default;

    ThreadData(vector<TUPLE_S> *const s, volatile bool *const dead, vector<int> *const skyline,
               const uint32_t thread_id, const uint32_t num_of_threads): S(s),
                                                                         dead(dead),
                                                                         skyline(skyline),
                                                                         thread_id(thread_id),
                                                                         num_of_threads(num_of_threads) {
    }
};

BSkyTree::BSkyTree(const uint32_t num_threads,
                   const uint32_t n, const uint32_t d) : num_threads_(num_threads), n_(n), d_(d) {
    skyline_.reserve(1024);
    eqm_.reserve(1024);
}

void BSkyTree::Init(float **dataset) {
    data_.reserve(n_);
    for (uint32_t i = 0; i < n_; i++) {
        TUPLE t;
        t.pid = i;
        memcpy(t.elems, dataset[i], sizeof(float) * NUM_DIMS);
        data_.push_back(TUPLE_S(t, -1));
    }
}

BSkyTree::~BSkyTree() {
    skyline_.clear();
    data_.clear();
    eqm_.clear();
}

vector<int> BSkyTree::Execute() {
    //  initProfiler();
    BSkyTreeS_ALGO();

    //  printProfile();

    // Add missing points from "equivalence matrix"
    skyline_.insert(skyline_.end(), eqm_.begin(), eqm_.end());

    return skyline_;
}

void BSkyTree::BSkyTreeS_ALGO() {
    //  initProfiler();
    SelectBalanced(); // pivot selection in the data_
    //  updateProfiler( "pivot" );

    DoPartitioning(); // mapping points to binary vectors representing sub-regions
    //  updateProfiler( "partitioning" );

    //  initProfiler();
    vector<TUPLE_S> &S = data_; // Alias
    volatile bool *dead = new bool[S.size()];
    memset(const_cast<bool *>(dead), 0, sizeof(bool) * S.size());

    head = 1;
    tail = S.size() - 1; // always points to the last (active) tuple

    // Initialize barrier
    pthread_barrier_init(&barrier, nullptr, num_threads_);

    // Create threads
    const auto tid = new pthread_t[num_threads_];
    const auto thread_data = new ThreadData [num_threads_]; // ThreadData *thread_data = new ThreadData [num_threads_];

    for (int i = 0; i < num_threads_; i++) {
        thread_data[i] = ThreadData(&S, dead, &skyline_, i, num_threads_);
        pthread_create(&tid[i], nullptr, parallel_execution, &thread_data[i]);
    }

    // Main thread waits for the termination of worker threads
    for (int i = 0; i < num_threads_; i++) {
        pthread_join(tid[i], nullptr);
    }

    // Destroy barrier
    pthread_barrier_destroy(&barrier);

    delete[] tid;
    delete[] thread_data;
    delete[] dead;
}

/*
 * Partitions the data using the pivot point (data_[0]) by
 * assigning partition bitmap to each tuple. Also, removes
 * the points that are pruned (ALL_ONES partition).
 */
void BSkyTree::DoPartitioning() {
    const uint32_t pruned = SHIFTS[NUM_DIMS] - 1;
    const TUPLE_S &pivot = data_[0];
    for (uint32_t i = 1; i < data_.size(); ++i) {
        if (EqualityTest(pivot, data_[i])) {
            eqm_.push_back(data_[i].pid);
            data_[i] = data_.back();
            data_.pop_back();
            continue;
        }
        const uint32_t lattice = DT_bitmap_dvc(data_[i], pivot);
        if (lattice < pruned) {
            assert(!DominateLeft( pivot, data_[i] ));
            data_[i].partition = lattice;
        } else {
            data_[i] = data_.back();
            data_.pop_back();
        }
    }
}

/*
 * Chooses a pivot based on minimum range. The chosen pivot
 * is a skyline point.
 *
 * In addition to that, removes points from data_ that are
 * dominated by the (current) pivot point.
 *
 * The pivot point is stored in data_[0].
 */
void BSkyTree::SelectBalanced() {
    const vector<float> min_list(NUM_DIMS, 0.0);
    const vector<float> max_list(NUM_DIMS, 1.0);

    constexpr uint32_t head = 0;
    uint32_t tail = data_.size() - 1, cur_pos = 1;
    const float *hvalue = data_[head].elems;

    const vector<float> range_list = SetRangeList(min_list, max_list);
    float min_dist = ComputeDistance(hvalue, min_list, range_list);

    while (cur_pos <= tail) {
        const float *cvalue = data_[cur_pos].elems;

        const uint32_t dtest = DominanceTest(data_[head], data_[cur_pos]);
        if (dtest == DOM_LEFT) {
            data_[cur_pos] = data_[tail];
            data_.pop_back();
            tail--;
        } else if (dtest == DOM_RIGHT) {
            data_[head] = data_[cur_pos];
            data_[cur_pos] = data_[tail];
            data_.pop_back();
            tail--;

            hvalue = data_[head].elems;
            min_dist = ComputeDistance(hvalue, min_list, range_list);
            cur_pos = 1; // THIS IS THE SAME BUG AS IN QSkyCube: cur_pos is not reseted
        } else {
            assert(dtest == DOM_INCOMP);
            const float cur_dist = ComputeDistance(cvalue, min_list, range_list);

            if (cur_dist < min_dist) {
                if (EvaluatePoint(cur_pos)) {
                    std::swap(data_[head], data_[cur_pos]);

                    // hvalue = data_[head].elems;
                    min_dist = cur_dist;
                    cur_pos++;
                } else {
                    data_[cur_pos] = data_[tail];
                    data_.pop_back();
                    tail--;
                }
            } else
                cur_pos++;
        }
    }
}

vector<float> BSkyTree::SetRangeList(const vector<float> &min_list,
                                     const vector<float> &max_list) {
    vector<float> range_list(NUM_DIMS, 0);
    for (uint32_t d = 0; d < NUM_DIMS; d++)
        range_list[d] = max_list[d] - min_list[d];

    return range_list;
}

/*
 * Note that here we do not need to do normalization (we assume
 * the data is pre-normalized). Though, we do (TODO!)
 */
float BSkyTree::ComputeDistance(const float *value,
                                const vector<float> &min_list, const vector<float> &range_list) {
    float min_d;

    float max_d = min_d = (value[0] - min_list[0]) / range_list[0];
    for (uint32_t d = 1; d < NUM_DIMS; d++) {
        const float norm_value = (value[d] - min_list[d]) / range_list[d];
        if (min_d > norm_value)
            min_d = norm_value;
        else if (max_d < norm_value)
            max_d = norm_value;
    }

    return max_d - min_d;
}

/*
 * Checks if the point dataset[pos] is not dominated by any of points
 * before pos (dataset[0..pos-1]).
 *
 * Note that here we can remove additionally dominated points, but the
 * code does not do it (the paper suggests, though).
 */
bool BSkyTree::EvaluatePoint(const uint32_t pos) const {
    const TUPLE &cur_tuple = data_[pos];
    for (uint32_t i = 0; i < pos; ++i) {
        const TUPLE &prev_value = data_[i];
        if (DominatedLeft(cur_tuple, prev_value))
            return false;
    }

    return true;
}

void *BSkyTree::parallel_execution(void *arg) {
    const auto *thread_data = static_cast<ThreadData *>(arg); // ThreadData *thread_data2 = (ThreadData *) arg;
    const auto thread_id = thread_data->thread_id;
    const auto num_of_threads = thread_data->num_of_threads;
    const auto S = thread_data->S;
    const auto dead = thread_data->dead;

    while (head < tail) {
        uint32_t htail = head + BSKYTREE_ACCUM - 1 < tail ? head + BSKYTREE_ACCUM - 1 : tail;
        const uint32_t p_head = head, p_tail = tail, p_htail = htail;
        uint32_t const num_heads = p_htail - p_head + 1;
        uint32_t const chunk = (num_heads + num_of_threads - 1) / num_of_threads;
        uint32_t const start = p_head + thread_id * chunk;
        uint32_t end = p_head + (thread_id + 1) * chunk - 1;
        end = end <= p_htail ? end : p_htail;

        for (uint32_t th = start; th <= end; th++) {
            // th -> temporal head
            uint32_t cur = p_htail + 1;
            while (cur <= p_tail) {
                if (dead[cur]) {
                    ++cur;
                    continue;
                }

                // Check for partial dominance:
                if (((*S)[th].partition & (*S)[cur].partition) == (*S)[th].partition
                    || ((*S)[th].partition & (*S)[cur].partition) == (*S)[cur].partition) {
                    const int dt_test = DT_dvc((*S)[th], (*S)[cur]);
                    if (dt_test == DOM_LEFT) {
                        dead[cur] = true;
                        //          S[cur++] = S[tail--]; no compression because of multi-threading
                        cur++;
                    } else if (dt_test == DOM_RIGHT) {
                        dead[cur] = true;
                        (*S)[th] = (*S)[cur];
                        cur = p_htail + 1;
                        //          S[cur] = S[tail--]; no compression because of multi-threading
                    } else {
                        cur++; // Point-level incomparability
                    }
                } else {
                    cur++; // Region-level incomparability
                }
            } // with cur we did one pass (till tail)
        } //  (all temporal heads were processed)

        pthread_barrier_wait(&barrier);

        if (thread_id == 0) {
            for (uint32_t th = head; th <= htail; ++th) {
                // th -> temporal head
                uint32_t c = th + 1;
                while (c <= htail) {
                    if ((*S)[th].pid == (*S)[c].pid) {
                        dead[htail] = true;
                        (*S)[c] = (*S)[htail--];
                    } else {
                        const uint32_t dt_test = DT_dvc((*S)[th], (*S)[c]);
                        if (dt_test == DOM_LEFT) {
                            dead[htail] = true;
                            (*S)[c] = (*S)[htail--];
                        } else if (dt_test == DOM_RIGHT) {
                            (*S)[th] = (*S)[c];
                            dead[htail] = true;
                            (*S)[c] = (*S)[htail--];
                            c = th + 1;
                        } else {
                            c++; // two heads are incomparable
                        }
                    }
                }
            }

            // heads are processed, so update it:
            head = htail + 1;
            //    updateProfiler( "sequential head processing" );

            // Compress by removing dead tuples:
            uint32_t head_dead = head; // will store 1st dead tuple from head & upwards
            uint32_t tail_alive = tail; // will store 1st alive tuple from tail & downwards
            while (head_dead < tail_alive) {
                while (!dead[head_dead] && head_dead < tail_alive) ++head_dead;

                while (dead[tail_alive]) --tail_alive;

                if (tail_alive > head_dead) {
                    dead[head_dead] = false;
                    dead[tail_alive] = true;
                    (*S)[head_dead++] = (*S)[tail_alive--];
                }
            }
            // Update tail after compression:
            tail = tail_alive;
            while (dead[tail]) --tail;
            //    updateProfiler( "sequential compression" );
        }

        pthread_barrier_wait(&barrier);
    } // Skyline computed!

    if (thread_id == 0) for (uint32_t i = 0; i <= tail; ++i) thread_data->skyline->push_back((*S)[i].pid);
    return nullptr;
}
