#include <bits/stdc++.h>
#include <sp_tree.hpp>

typedef uint32_t K;

const K N = (1<<20), Q = (1<< 18); // <- change these
SP_Tree<K> sp_tree;
int main() {
    // N = 123456;
    // printf("N = %d, Q = %d\n", N, Q);

    std::mt19937 rng(0);
    //
    std::vector<K> a;
    std::vector<K> q;
    //
    a.resize(N);
    q.resize(Q);

    for (K i = 0; i < N; i++) {
        a[i] = i * 2;
        // std::cerr << a[i] << " ";
    }
    std::cerr << std::endl;

    std::sort(a.begin(), a.end());
    sp_tree.build(a);
    // sp_tree.save_tree("./sp_tree.bin");

    for (K i = 0; i < Q; i++) {
        q[i] = i * 2 + 1;
        std::cerr << q[i] << " " << sp_tree.lower_bound_index(q[i]) << std::endl;

    }
    // sp_tree.load_tree("./sp_tree.bin");

    return 0;
}