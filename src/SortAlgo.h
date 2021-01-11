/******************************************************************************
 * src/SortAlgo.h
 *
 * Implementations of many sorting algorithms.
 *
 * Note that these implementations may not be as good/fast as possible. Some
 * are modified so that the visualization is more instructive.
 *
 * Futhermore, some algorithms are annotated using the mark() and watch()
 * functions from SortArray. These functions add colors to the illustratation
 * and thereby makes the algorithm's visualization easier to explain.
 *
 ******************************************************************************
 * Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef SORTALGO_H
#define SORTALGO_H

#include <wx/string.h>
#include "SortArray.h"

// *** List of Sorting Algorithms

struct AlgoEntry
{
    wxString name;
    void (*func)(class SortArray&);
    // maximum item count for test runs
    unsigned int max_testsize;
    // count inversions if n <= limit
    unsigned int inversion_count_limit;
    wxString text;
};

extern const struct AlgoEntry g_algolist[];
extern const size_t g_algolist_size;
extern const struct AlgoEntry* g_algolist_end;

// *** Sorting Algorithms

void SelectionSort(class SortArray& a);
void InsertionSort(class SortArray& a);
void BinaryInsertionSort(class SortArray& a);

void MergeSort(class SortArray& a);
void MergeSortIterative(class SortArray& a);

wxArrayString QuickSortPivotText();

enum QuickSortPivotType { PIVOT_FIRST, PIVOT_LAST, PIVOT_MID, PIVOT_RANDOM, PIVOT_MEDIAN3, PIVOT_MOM};
extern QuickSortPivotType g_quicksort_pivot;

void QuickSortLR(class SortArray& a);
void QuickSortLL(class SortArray& a);
void QuickSortTernaryLR(class SortArray& a);
void QuickSortTernaryLL(class SortArray& a);
void QuickSortDualPivot(class SortArray& a);

void BubbleSort(class SortArray& a);
void CocktailShakerSort(class SortArray& a);
void CombSort(class SortArray& a);
void GnomeSort(class SortArray& a);
void OddEvenSort(class SortArray& a);

void ShellSort(SortArray& a);
void HeapSort(class SortArray& a);
void SmoothSort(class SortArray& a);

void BitonicSortRe(SortArray& a);
void BitonicSortIt(SortArray& a);
void BatcherSortNetworkRe(SortArray& a);
void BatcherSortNetworkIt(SortArray& a);

void RadixSortLSD(class SortArray& a);
void RadixSortMSD(class SortArray& a);

void StlSort(class SortArray& a);
void StlStableSort(class SortArray& a);
void StlHeapSort(class SortArray& a);

void TimSort(class SortArray& a);
void WikiSort(class SortArray& a);

void BogoSort(class SortArray& a);
void BozoSort(class SortArray& a);
void StoogeSort(class SortArray& a);

void CycleSort(class SortArray& a);
void SlowSort(class SortArray& a);

// ****************

void CircleSort(class SortArray& a);
void IterativeCircleSort(class SortArray& a);
void CircloidSort(class SortArray& a);

void KiyomiSort(class SortArray& a);
void ChinottoSort(class SortArray& a);
void KinnowSort(class SortArray& a);
void MandarinOrangeSort(class SortArray& a);

void BubbleScanQuicksort(class SortArray& a);
void StableQuickSort(class SortArray& a);
void OptiStableQuickSort(class SortArray& a);
void QuickSortLLL(class SortArray& a);
void IterativeQuickSortLL(class SortArray& a);
void IterativeQuickSortLR(class SortArray& a);
void IterativeLinkedQuickSortLL(class SortArray& a);
void IterativeLinkedQuickSortLR(class SortArray& a);

void SmarterBubbleSort(class SortArray& a);

void ModuloSort(class SortArray& a);
void RotateRadixLSD(class SortArray& a);
void RotateRadixMSD(class SortArray& a);
void ProxmapSort(class SortArray& a);

void FloatSort(class SortArray& a);
void OddEvenBase3(class SortArray& a);

void StableSelectionSort(class SortArray& a);
void ShoveSort(class SortArray& a);

void PairwiseSortingIt(class SortArray& a);
void PairwiseSortingRe(class SortArray& a);
void BoseNelsonSortingIt(class SortArray& a);
void BoseNelsonSortingRe(class SortArray& a);
void BalancedSortingNetwork(class SortArray& a);

void IterativeStoogeSort(class SortArray& a);
void QuadStoogeSort(class SortArray& a);
void StableStoogeSort(class SortArray& a);
void HyperStoogeSort(class SortArray& a);
void RadixLSDStoogeSort(class SortArray& a);
void RadixMSDStoogeSort(class SortArray& a);
void YSlowSort(class SortArray& a);
void GravitySort(class SortArray& a);
void LessBogoSort(class SortArray& a);
void AwkwardSort(class SortArray& a);
void SnuffleSort(class SortArray& a);

void GrateSort(class SortArray& a);
void ReverseGrateSort(class SortArray& a);
void CocktailGrateSort(class SortArray& a);
void RoomSort(class SortArray& a);
void SlopeSort(class SortArray& a);
void BalanceSort(class SortArray& a);
void WiggleSort(class SortArray& a);
void BogoBogoBogoSort(class SortArray& a);
void SlideSort(class SortArray& a);
void BinarySort(class SortArray& a);

void JumpDownSort(class SortArray& a);
void CocktailShellSort(class SortArray& a);
void TernarySlowSort(class SortArray& a);
void OddEvenBogoSort(class SortArray& a);

void PancakeSort(class SortArray& a);
void QuasiPancakeSort(class SortArray& a);
void BinaryQuasiPancakeSort(class SortArray& a);

// void ImprovedWeaveMergeSort(class SortArray& a);
void BingoSort(class SortArray& a);
void ReverseInsertionSort(class SortArray& a);
void OddEvenTransMergeSort(class SortArray& a);
void SpaghettiSort(class SortArray& a);

void UnbalancedTreeSort(class SortArray& a);
void FlippedMinHeapSort(class SortArray& a);
void CheckerboardHeapSort(class SortArray& a);
void TriangularHeapSort(class SortArray& a);

// ****************

void InsertionSortExtra(class SortArray& a, int lo, int hi);
ssize_t QuickSortSelectPivot(class SortArray& a, ssize_t lo, ssize_t hi);
bool CheckSorted(class SortArray& a);

// ****************************************************************************
// *** Iterator Adapter

// iterator based on http://zotu.blogspot.de/2010/01/creating-random-access-iterator.html

class MyIterator : public std::iterator<std::random_access_iterator_tag, ArrayItem>
{
protected:
    SortArray*  m_array;
    size_t      m_pos;

public:
    typedef std::iterator<std::random_access_iterator_tag, ArrayItem> base_type;

    typedef std::random_access_iterator_tag iterator_category;

    typedef base_type::value_type value_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::reference reference;
    typedef base_type::pointer pointer;

    MyIterator() : m_array(NULL), m_pos(0) {}

    MyIterator(SortArray* A, size_t p) : m_array(A), m_pos(p) {}

    MyIterator(const MyIterator& r) : m_array(r.m_array), m_pos(r.m_pos) {}

    MyIterator& operator=(const MyIterator& r)
    { m_array = r.m_array, m_pos = r.m_pos; return *this; }

    MyIterator& operator++()
    { ++m_pos; return *this; }

    MyIterator& operator--()
    { --m_pos; return *this; }

    MyIterator operator++(int)
    { return MyIterator(m_array, m_pos++); }

    MyIterator operator--(int)
    { return MyIterator(m_array, m_pos--); }

    MyIterator operator+(const difference_type& n) const
    { return MyIterator(m_array, m_pos + n); }

    MyIterator& operator+=(const difference_type& n)
    { m_pos += n; return *this; }

    MyIterator operator-(const difference_type& n) const
    { return MyIterator(m_array, m_pos - n); }

    MyIterator& operator-=(const difference_type& n)
    { m_pos -= n; return *this; }

    reference operator*() const
    { return m_array->get_mutable(m_pos); }

    pointer operator->() const
    { return &(m_array->get_mutable(m_pos)); }

    reference operator[](const difference_type& n) const
    { return m_array->get_mutable(m_pos + n); }

    bool operator==(const MyIterator& r)
    { return (m_array == r.m_array) && (m_pos == r.m_pos); }

    bool operator!=(const MyIterator& r)
    { return (m_array != r.m_array) || (m_pos != r.m_pos); }

    bool operator<(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos < r.m_pos) : (m_array < r.m_array)); }

    bool operator>(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos > r.m_pos) : (m_array > r.m_array)); }

    bool operator<=(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos <= r.m_pos) : (m_array <= r.m_array)); }

    bool operator>=(const MyIterator& r)
    { return (m_array == r.m_array ? (m_pos >= r.m_pos) : (m_array >= r.m_array)); }

    difference_type operator+(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos + r2.m_pos); }

    difference_type operator-(const MyIterator& r2) const
    { ASSERT(m_array == r2.m_array); return (m_pos - r2.m_pos); }
};


#endif // SORTALGO_H
