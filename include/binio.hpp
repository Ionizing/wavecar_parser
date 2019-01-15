#pragma once

#ifndef BINIO_H
#define BINIO_H

// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACK

#include <base.hpp>


// all the matrices or vectors are RowMajor
namespace ionizing {
  class BinIO {
  public:
    BinIO(const char* FileName);
    ~BinIO();
                            // Prevent 'copying' operations.
    BinIO(const BinIO&)            = delete;
    BinIO& operator=(const BinIO&) = delete;
                            // change pointer focus
    void seek(const int n);
    int getFileSize() const;

  template <typename T>
    T readElement();

  template <typename T>
    ColVecT<T> readVectorCol(const long size);

  template <typename T>
    RowVecT<T> readVectorRow(const long size);

  template <typename T>
    MatT<T> readMatrix(const long nRow, const long nCol);

  private:
    std::ifstream ifs;
    std::ofstream ofs;
    int _fileSize;

  private:
    void _openFile(const char* FileName);

  }; // end of class BinIO
};


namespace ionizing {

  template <typename T>
    T BinIO::readElement() {
      long len = sizeof(T);
      T x;
      ifs.read((char*) &x, len);
      return x;
    }

  template <typename T>
    MatT<T> BinIO::readMatrix(const long nRow, const long nCol) {
      Matrix<T, Dynamic, Dynamic, RowMajor> matrix(nRow, nCol);
      ifs.read((char*) matrix.data(), nRow * nCol * sizeof(T));
      return matrix;
    }

  template <typename T> 
    ColVecT<T> BinIO::readVectorCol(const long size) {
      return readMatrix<T>(size, 1);
    }
 
  template <typename T>
    RowVecT<T> BinIO::readVectorRow(const long size) {
      return readMatrix<T>(1, size);
    } 


}
#endif  // BINIO_H
