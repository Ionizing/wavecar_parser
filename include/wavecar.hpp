#pragma once

#ifndef WAVECAR_H
#define WAVECAR_H

#include <base.hpp>
#include <binio.hpp>

namespace ionizing {
  // type aliasing declared in base.hpp, statements here are helpful
  // to understand this code a year later.

  // using namespace Eigen;
  // using Vecd   = VectorXd;
  // using Veccd  = VectorXcd;
  // using Matd   = Matrix<double,               Dynamic, Dynamic>;
  // using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic>;
  // using Cubd   = Tensor<double,                              3, RowMajor>;
  // using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
  // using Mat33d = Matrix<double,                     3,       3>;
  // using MatX3d = Matrix<double,               Dynamic,       3>;


class WAVECAR {
public:
  // type defs
  struct Info {
    int  _recordLength,
         _nSpin,
         _precisionTag;
    bool _isDoubleType;
    int  _fileSize;
  };

  struct Header {
    int    _nKpoints,
           _nBands;
    double _enCut;
    Mat33d _latticeVectors,
           _reciprocalVectors;
    double _omega;
    double _eFermi;
  };




  WAVECAR(const char* FileName);
  WAVECAR(const WAVECAR&) = delete;
  WAVECAR operator=(const WAVECAR&) = delete;


  bool getPrecisionType() const;
                // returns recl, nspin, prectag, isdoubletype, filesize
  const Info&                  getInfo()              const;
                // returns header: nkpts, nbands, encut, Acell
  const Header&                getHeader()            const;
  const MatX3d&                getKVectors()          const;
  const ColVecT<Matcd>&        getBands()             const;
  const ColVecT<Matd>&         getFermiWeights()      const;
  const ColVecT<MatT<Veccd>>&  getComplexWaves()      const;
  const Mat33d&                getReciprocalVectors() const;
  const Mat33d&                getLatticeVectors()    const;
  const VectorXi&              getNPlaneWaves()       const;


  void printInfo(std::ostream& os) const;
  const bool checkIndex(const int  ispin,
                        const int  ikpoint, 
                        const int  iband) const;
  const Veccd getBandCoeff(const int   ispin,
                           const int   ikpoint,
                           const int   iband,
                           const bool  is_normed) const;
  
 

private:
  BinIO io;                           // binary io operators, non-copyable
  int   _iRec;                        // Record index corresponding to IREC in fileio.F

    Info   _info; 
    Header _header;

  // struct _body {
    VectorXi       _nPlaneWaves;      // int array, number of plane waves for each kpoint
    MatX3d         _kVectors;         // X * 3 array, k-vectors for each kpoint
    ColVecT<Matcd> _bands;            // nspin * nkpts * nbands
    ColVecT<Matd>  _fermiWeights;
    int            _maxOfNPlaneWaves;
                                      // _complexWaves 4 dimensions are: 
                                      //              (ispin, kpoint, nband, nplanewaves)
    ColVecT<MatT<Veccd>> _complexWaves;
  // };

private:
  void read_info();
  void read_header();
  void read_band();

}; // class WAVECAR


} // namespace ionizing

#endif // WAVECAR_H
