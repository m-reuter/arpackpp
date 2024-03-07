
#ifndef ARSPMAT_H
#define ARSPMAT_H

#include <cstring>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arerror.h"

/**
 * @brief A sparse matrix represented in compressed column storage.
 */
template<typename ARTYPE>
class ARSparseMatrix : public ARMatrix<ARTYPE> {

protected:

    int*    Ap;    // column pointers
    int*    Ai;    // row indices
    ARTYPE* Ax;    // values
    int*    Di;    // diagonal indices (allocated as needed)
    int     nzmax; // currently allocated space for non-zeros in Ai and Ax
    char    uplo;  // symmetric matrix given as upper 'U' or lower 'L'
    bool    owner; // if true, Ap, Ai and Ax were allocated by this class.

    void ClearMem();

    /**
     * @brief Expands a lower triangular sparse matrix.
     *
     * @param B The target matrix B
     * @param work Work array of size n+1
     *
     * @return Returns 0 on success. Otherwise the number of elements exceeds
     *      nzmax and the minimum required value for nzmax is returned.
     */
    int ExpandL(ARSparseMatrix<ARTYPE>& B, int* work);

    /**
     * @brief Expands an upper triangular sparse matrix.
     *
     * @param B The target matrix B
     * @param work Work array of size n+1
     *
     * @return Returns 0 on success. Otherwise the number of elements exceeds
     *      nzmax and the minimum required value for nzmax is returned.
     */
    int ExpandU(ARSparseMatrix<ARTYPE>& B, int* work);

public:

    /**
     * @brief Returns the column pointers array (size ncols()).
     */
    int* pcol() { return Ap; }

    /**
     * @brief Returns the row indices array (size nzeros()).
     */
    int* irow() { return Ai; }

    /**
     * @brief Returns the values array (size nzeros()).
     */
    ARTYPE* values() { return Ax; }

    /**
     * @brief Returns the number of non-zeros.
     */
    int nzeros() const { return Ap[this->n]; }

    /**
     * @brief Returns the maximum number of non-zeros.
     */
    int size() const { return nzmax; }

    bool IsUpper() const { return uplo == 'U'; }
    bool IsLower() const { return uplo == 'L'; }
    bool IsTriangular() const { return uplo == 'L' || uplo == 'U'; }

    /**
     * @brief Check whether the compressed column storage is valid.
     */
    bool Check() const;

    /**
     * @brief Copy the contents of given source matrix to this instance.
     *
     * @param other The source matrix
     */
    virtual void Copy(const ARSparseMatrix<ARTYPE>& other);

    /**
     * @brief Expands a lower or upper triangular sparse matrix.
     *
     * @param B The target matrix B
     */
    int Expand(ARSparseMatrix<ARTYPE>& B);

    /**
     * @brief Computes the positions of the diagonal elements of a sparse matrix
     *        and stores the indices in the Di array.
     *
     * @param update Force update of diagonal indices, even if already computed
     *
     * @return On first call, the number of diagonal elements present in the sparse
     *         matrix is returned. Subsequent calls will return -1.
     *
     * @remarks The value Di[i] is the index of the diagonal element A(i,i) in
     *          the arrays Ax and Ap of the sparse matrix. If no diagonal element
     *          is found, the entry is set to a negative int.
     */
    int DiagIndices(bool update = false);

    /**
     * @brief Adds a scalar to the diagonal entries of a sparse matrix A = A + s I
     *
     * @param value Scalar to add to the diagonal entries.
     *
     * @return The number of non-zero entries in the sparse matrix.
     *
     * @remarks The matrix A may be expanded to allow for additions of nonzero
     *    elements to previously non-existing diagonals. There is no checking
     *    as to whether there is enough space appended to the arrays Ax and Ap.
     *    If not sure, allow for n additional elements.
     */
    int AddDiag(ARTYPE value);

    /**
     * @brief Performs the operation C = A + s B.
     *
     * @param s Scalar factor for B
     * @param B The matrix to add
     * @param C The target matrix
     *
     * @return Returns 0 if ok, otherwise the number of elements in C exceeds nzmax
     */
    int Add(const ARTYPE s, const ARSparseMatrix<ARTYPE>& B, ARSparseMatrix<ARTYPE>& C);

    /**
     * @brief Performs the operation C = A + s B (special overload for real-valued, non-symmetric
     *        problems with complex shift).
     *
     * @param sr Real part of complex scalar factor for B
     * @param si Imaginary part of complex scalar factor for B
     * @param B  The real-valued matrix to add
     * @param C  The complex-valued target matrix
     *
     * @return Returns 0 on success, otherwise the number of elements in C exceeds nzmax
     */
    int Add(const ARTYPE sr, const ARTYPE si, const ARSparseMatrix<ARTYPE>& B, ARSparseMatrix<arcomplex<ARTYPE>>& C);

    /**
     * @brief Gets the exact number of nonzero elements in A + B.
     *
     * @param B The matrix to add
     * @param count Array of length ncol containing the degrees (i.e. the
     *              number of non-zeros in each column of the matrix A + B
     * @param work Work array of length nrow
     *
     * @return Total number of nonzero elements in A + B
     */
    int PrepareAdd(const ARSparseMatrix<ARTYPE>& B, int* count, int* work);

    void MultMv(ARTYPE* v, ARTYPE* w);


    ARSparseMatrix(int nrows, int ncols, int nz, char uplo = '*')
        : ARMatrix<ARTYPE>(nrows, ncols), nzmax(nz), uplo(uplo), owner(true)
    {
        Ap = new int[ncols + 1];
        Ai = new int[nz];
        Ax = new ARTYPE[nz];

        // Make sure the CSC arrays are initialized correctly.
        std::fill(Ap, Ap + ncols + 1, 0);
        std::fill(Ai, Ai + nz, 0);

        Di = nullptr;
    }

    ARSparseMatrix(int nrows, int ncols, int* &ap, int* &ai, ARTYPE* &ax, int nzmax = -1, char uplo = '*', bool owner = false)
        : ARMatrix<ARTYPE>(nrows, ncols), Ap(ap), Ai(ai), Ax(ax), uplo(uplo), owner(owner)
    {
      if (ap == nullptr)
      {
        // throw?
      }

      this->nzmax = (nzmax < 0) ? Ap[ncols] : nzmax;

      Di = nullptr;
    }

    ARSparseMatrix(const ARSparseMatrix<ARTYPE>& other) { Copy(other); }
    // Copy constructor.

    virtual ~ARSparseMatrix() { ClearMem(); }
    // Destructor.

    ARSparseMatrix& operator=(const ARSparseMatrix<ARTYPE>& other);
    // Assignment operator.

};

template<typename ARTYPE>
inline void ARSparseMatrix<ARTYPE>::ClearMem()
{
    if (owner)
    {
        if (Ap) { delete[] Ap; Ap = nullptr; }
        if (Ai) { delete[] Ai; Ai = nullptr; }
        if (Ax) { delete[] Ax; Ax = nullptr; }
    }

    if (Di) { delete[] Di; Di = nullptr; }
}

template<typename ARTYPE>
bool ARSparseMatrix<ARTYPE>::Check() const
{
    int i, j, k;
    
    int m = this->m;
    int n = this->n;

    // Checking if column pointers are in ascending order.

    i = 0;
    while (i != n && Ap[i] <= Ap[i + 1]) i++;

    if (i != n) return false;

    // Checking if row indices are in order and within bounds.

    for (i = 0; i != n; i++)
    {
        j = Ap[i];
        k = Ap[i + 1] - 1;

        if (j <= k)
        {
            if (uplo == 'L')
            {
                if (Ai[j] < i || Ai[k] >= m) return false;
            }
            else if (uplo == 'U')
            {
                if (Ai[j] < 0 || Ai[k] > i) return false;
            }
            else
            {
                if (Ai[j] < 0 || Ai[k] >= m) return false;
            }

            while ((j != k) && (Ai[j] < Ai[j + 1])) j++;

            if (j != k) return false;
        }
    }

    return true;
}

template <typename ARTYPE>
inline void ARSparseMatrix<ARTYPE>::Copy(const ARSparseMatrix<ARTYPE>& other)
{
    int n = this->n;

    if (n != other.n)
    {
        throw ArpackError(ArpackError::INCOMPATIBLE_SIZES, "ARSparseMatrix::Copy");
    }

    int nnz = other.nzeros();

    if (nzmax < nnz)
    {
        throw ArpackError(ArpackError::INCOMPATIBLE_SIZES, "ARSparseMatrix::Copy");
    }

    uplo = other.uplo;

    std::memcpy(Ap, other.Ap, (n + 1) * sizeof(int));
    std::memcpy(Ai, other.Ai, nnz * sizeof(int));
    std::memcpy(Ax, other.Ax, nnz * sizeof(ARTYPE));

    if (other.Di)
    {
        if (!Di)
        {
            Di = new int[n];
        }
        std::memcpy(Di, other.Di, n * sizeof(int));
    }
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::DiagIndices(bool update)
{
    // The diagonal indices have already been computed.
    if (Di && !update) return -1;

    int j, end;
    int count = 0;

    int n = this->n;

    if (!Di)
    {
        Di = new int[n];
    }

    for (int i = 0; i < n; i++)
    {
        j = Ap[i];
        end = Ap[i + 1];

        // Skip super-diagonals.
        while ((Ai[j] < i) && (j < end))
            j++;

        // Is diagonal entry?
        if ((Ai[j] == i) && (j < end))
        {
            Di[i] = j;
            count++;
        }
        else
        {
            Di[i] = -j - 1;
        }
    }

    return count;
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::AddDiag(ARTYPE value)
{
    int i, j, k, k0, start, end;
    bool test;
    int n = this->n;

    // Ensure diagonal indices are computed.
    DiagIndices();

    // Number of missing diagonals.
    int icount = 0;

    for (i = 0; i < n; i++)
    {
        if (Di[i] < 0)
        {
            ++icount;
        }
        else
        {
            Ax[Di[i]] += value;
        }
    }

    // Return if no diagonal elements are missing.
    if (icount == 0)
    {
        return Ap[n];
    }

    // Shift the nonzero elements if needed, to allow for created diagonal elements.

    k0 = Ap[n] + icount;

    // Copy columns backward.
    for (i = n - 1; i >= 0; --i)
    {
        // Go through column i.
        start = Ap[i];
        end = Ap[i + 1] - 1;
        Ap[i + 1] = k0;
        test = Di[i] < 0;
        for (k = end; k >= start; --k)
        {
            j = Ai[k];
            if (test && j < i)
            {
                test = false;
                k0--;
                Ax[k0] = value;
                Ai[k0] = i;
                Di[i] = k0;
            }
            k0--;
            Ax[k0] = Ax[k];
            Ai[k0] = j;
        }

        // Diagonal element has not been added yet.
        if (test)
        {
            k0--;
            Ax[k0] = value;
            Ai[k0] = i;
            Di[i] = k0;
        }
    }
    Ap[0] = k0;
    return Ap[n];
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::Add(const ARTYPE s, const ARSparseMatrix<ARTYPE> &B,
    ARSparseMatrix<ARTYPE>& C)
{
    if (s == (ARTYPE)0.0)
    {
        C.Copy(*this);
        return 0;
    }

    int i, rowa, rowb, kc, ka, kb, kamax, kbmax, nzmax;

    int nrow = this->m;
    int ncol = this->n;

    auto Bp = B.Ap;
    auto Bi = B.Ai;
    auto Bx = B.Ax;

    auto Cp = C.Ap;
    auto Ci = C.Ai;
    auto Cx = C.Ax;

    nzmax = C.nzmax;

    Cp[0] = kc = 0;

    for (i = 0; i < ncol; i++)
    {
        ka = Ap[i];
        kb = Bp[i];
        kamax = Ap[i + 1];
        kbmax = Bp[i + 1];

        while (ka < kamax || kb < kbmax)
        {
            rowa = (ka < kamax) ? Ai[ka] : nrow;
            rowb = (kb < kbmax) ? Bi[kb] : nrow;

            if (rowa == rowb)
            {
                Cx[kc] = Ax[ka] + s * Bx[kb];
                Ci[kc] = rowa;
                ka++;
                kb++;
                kc++;
            }
            else if (rowa < rowb)
            {
                Ci[kc] = rowa;
                Cx[kc] = Ax[ka];
                ka++;
                kc++;
            }
            else if (rowa > rowb)
            {
                Ci[kc] = rowb;
                Cx[kc] = s * Bx[kb];
                kb++;
                kc++;
            }

            if (kc > nzmax)
            {
                return i;
            }
        }
        Cp[i + 1] = kc;
    }
    return 0;
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::Add(const ARTYPE sr, const ARTYPE si,
    const ARSparseMatrix<ARTYPE>& B, ARSparseMatrix<arcomplex<ARTYPE>>& C)
{
    auto Cp = C.pcol();
    auto Ci = C.irow();
    auto Cx = C.values();

    int nrow = this->m;
    int ncol = this->n;

    if (sr == (ARTYPE)0.0 && si == (ARTYPE)0.0)
    {
        int nnz = nzeros();

        std::memcpy(Cp, Ap, (ncol + 1) * sizeof(int));
        std::memcpy(Ci, Ai, nnz * sizeof(int));
        for (int i = 0; i < nnz; i++) Cx[i] = Ax[i];

        return 0;
    }

    int i, rowa, rowb, kc, ka, kb, kamax, kbmax, nzmax;

    auto Bp = B.Ap;
    auto Bi = B.Ai;
    auto Bx = B.Ax;

    nzmax = C.size();

    arcomplex<ARTYPE> s = arcomplex<ARTYPE>(sr, si);

    Cp[0] = kc = 0;

    for (i = 0; i < ncol; i++)
    {
        ka = Ap[i];
        kb = Bp[i];
        kamax = Ap[i + 1];
        kbmax = Bp[i + 1];

        while (ka < kamax || kb < kbmax)
        {
            rowa = (ka < kamax) ? Ai[ka] : nrow;
            rowb = (kb < kbmax) ? Bi[kb] : nrow;

            if (rowa == rowb)
            {
                Cx[kc] = Ax[ka] + s * Bx[kb];
                Ci[kc] = rowa;
                ka++;
                kb++;
                kc++;
            }
            else if (rowa < rowb)
            {
                Ci[kc] = rowa;
                Cx[kc] = Ax[ka];
                ka++;
                kc++;
            }
            else if (rowa > rowb)
            {
                Ci[kc] = rowb;
                Cx[kc] = s * Bx[kb];
                kb++;
                kc++;
            }

            if (kc > nzmax)
            {
                return i;
            }
        }
        Cp[i + 1] = kc;
    }
    return 0;
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::PrepareAdd(const ARSparseMatrix<ARTYPE>& B, int *count, int *work)
{
    int i, j, k, nz, last, end;

    auto Bp = B.Ap;
    auto Bi = B.Ai;

    int nrow = this->m;
    int ncol = this->n;

    for (i = 0; i < nrow; i++)
        work[i] = 0;
    for (i = 0; i < ncol; i++)
        count[i] = 0;

    for (i = 0; i < ncol; i++)
    {
        nz = 0;

        last = -1;

        // Columns of A

        end = Ap[i + 1];
        for (j = Ap[i]; j < end; j++)
        {
            k = Ai[j];

            nz++;
            work[k] = last;
            last = k;
        }

        // Columns of B

        end = Bp[i + 1];
        for (j = Bp[i]; j < end; j++)
        {
            k = Bi[j];

            if (work[k] == 0)
            {
                nz++;
                work[k] = last;
                last = k;
            }
        }

        count[i] = nz;

        // Reset work array

        for (j = 0; j < nz; j++)
        {
            k = work[last];
            work[last] = 0;
            last = k;
        }
    }

    nz = 0;

    for (i = 0; i < ncol; i++)
    {
        nz += count[i];
    }

    return nz;
}

template <typename ARTYPE>
inline void ARSparseMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{
    int    i, j;
    ARTYPE t;

    // Quitting the function if A was not defined.

    if (Ap == nullptr) {
        throw ArpackError(ArpackError::DATA_UNDEFINED, "ARSparseMatrix::MultMv");
    }

    // Determining w = M.v.

    for (i = 0; i != this->m; i++) w[i] = (ARTYPE)0.0;

    for (i = 0; i != this->n; i++)
    {
        t = v[i];
        for (j = Ap[i]; j != Ap[i + 1]; j++)
        {
            w[Ai[j]] += t * Ax[j];
        }
    }

}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::Expand(ARSparseMatrix<ARTYPE>& B)
{
    int i = 0;

    int* work = new int[this->n + 1];

    if (uplo == 'U')
    {
        i = ExpandU(B, work);

    }
    else if (uplo == 'L')
    {
        i = ExpandL(B, work);
    }

    // S = symmetric (not used anywhere, but might be useful).
    B.uplo = 'S';

    delete[] work;

    return i;
}

/* Helper function to conjugate values in the expanded matrix. */
template <typename T> inline T spmat_conj(T val) { return val; }
template <> inline arcomplex<float> spmat_conj(arcomplex<float> val) { return std::conj(val); }
template <> inline arcomplex<double> spmat_conj(arcomplex<double> val) { return std::conj(val); }

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::ExpandL(ARSparseMatrix<ARTYPE>& B, int *work)
{
    int i, j, k, nnz, ipos, end, nzmax;

    int n = this->n;

    auto Bp = B.Ap;
    auto Bi = B.Ai;
    auto Bx = B.Ax;

    nzmax = B.nzmax;

    for (i = 0; i < n + 1; i++)
    {
        Bp[i] = 0;
    }

    // Compute the number of elements in each row of strict lower part.

    for (i = 0; i < n; i++)
    {
        end = Ap[i + 1];
        for (k = Ap[i]; k < end; k++)
        {
            j = Ai[k];
            if (j > i)
            {
                Bp[j + 1]++;
            }
        }
    }

    // Find positions of first elements of output matrix.
    for (i = 0; i < n; i++)
    {
        Bp[i + 1] = Bp[i] + Bp[i + 1] + (Ap[i + 1] - Ap[i]);
    }

    // Enough storage?
    nnz = Bp[n];
    if (nnz > nzmax)
    {
        return nnz;
    }

    // Copy lower part.

    for (i = 0; i < n; i++)
    {
        end = Ap[i];
        k = Bp[i + 1] - 1;
        for (j = Ap[i + 1] - 1; j >= end; j--)
        {
            Bx[k] = Ax[j];
            Bi[k] = Ai[j];
            k--;
        }
        work[i] = Bp[i];
    }
    work[n] = Bp[n];

    // Now copy upper part (backwards). Go through the structure of Bx, Bp, Bi
    // that has already been copied (lower part). work(i) is the position of the
    // next free location in column i for Bx, Bp.

    for (i = 0; i < n; i++)
    {
        // i-th column is now in Bx, Bp, Bi structure -- upper half part
        end = Bp[i];
        for (k = Bp[i + 1] - 1; k >= end; k--)
        {
            j = Bi[k];
            if (j <= i)
            {
                break;
            }
            ipos = work[j];
            Bx[ipos] = spmat_conj(Bx[k]);
            Bi[ipos] = i;
            work[j]++;
        }
    }

    return 0;
}

template <typename ARTYPE>
inline int ARSparseMatrix<ARTYPE>::ExpandU(ARSparseMatrix<ARTYPE>& B, int *work)
{
    int i, j, k, nnz, ipos, end, nzmax;

    int n = this->n;

    auto Bp = B.Ap;
    auto Bi = B.Ai;
    auto Bx = B.Ax;

    nzmax = B.nzmax;

    for (i = 0; i < n + 1; i++)
    {
        Bp[i] = 0;
    }

    // Compute the number of elements in each row of strict upper part.

    for (i = 0; i < n; i++)
    {
        end = Ap[i + 1];
        for (k = Ap[i]; k < end; k++)
        {
            j = Ai[k];
            if (j < i)
            {
                Bp[j + 1]++;
            }
        }
    }

    // Find positions of first elements of output matrix.
    for (i = 0; i < n; i++)
    {
        Bp[i + 1] = Bp[i] + Bp[i + 1] + (Ap[i + 1] - Ap[i]);
    }
    // Enough storage?
    nnz = Bp[n];
    if (nnz > nzmax)
    {
        return nnz;
    }

    // Copy upper part (backwards).

    for (i = n - 1; i >= 0; --i)
    {
        end = Ap[i + 1];
        k = Bp[i];
        for (j = Ap[i]; j < end; j++)
        {
            Bx[k] = Ax[j];
            Bi[k] = Ai[j];
            k++;
        }
        work[i] = k;
    }
    work[n] = Bp[n];

    // Now copy lower part. Go through the structure of Bx, Bp, Bi that has
    // already been copied (lower part). work(i) is the position of the next
    // free location in column i for Bx, Bp.

    for (i = 0; i < n; i++)
    {
        // i-th column is now in Bx, Bp, Bi structure -- lower half part
        end = Bp[i + 1];
        for (k = Bp[i]; k < end; k++)
        {
            j = Bi[k];
            if (j >= i)
            {
                break;
            }
            ipos = work[j];
            Bx[ipos] = spmat_conj(Bx[k]);
            Bi[ipos] = i;
            work[j]++;
        }
    }

    return 0;
}

#endif
