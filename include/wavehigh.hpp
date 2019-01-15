#pragma once
#ifndef WAVEHIGH_H
#define WAVEHIGH_H

#include <wavecar.hpp>
#include <constants.hpp>
#include <fft.hpp>
#include <stringops.hpp>

namespace ionizing {

class WAVEHIGH : public WAVECAR {
  public:
    WAVEHIGH(const char* wavecar);
    WAVEHIGH(const WAVEHIGH&)           = delete;
    WAVEHIGH operator=(const WAVEHIGH&) = delete;
    void            plotBand()                          const;
    const Vecd&     getKPath()                          const;
    const Vector3i& getNGrid()                          const;
    const MatrixX3i getGVectors(const int ikpoint)      const;

    const Cubcd getKSWave(const int       ispin,
                          const int       ikpoint,
                          const int       iband,
                                MatrixX3i gvectors,
                                Vector3i  ngrid,
                          const Veccd&    coeff_vec,
                          const int       rescale,
                          const bool      is_norm) const;

    void saveAsVesta(const Cubcd& phi,
                     const char*  POSCAR,
                     const char*  prefix,
                     const bool   is_real) const;

    void plotWave(const int ispin,
                  const int ikpoint,
                  const int iband,
                  const char* poscar,
                  const char* prefix,
                  const bool is_real) const;





    // void plotDos() const;    // TODO
  private:
    char       _gammaHalf;
    bool       _isGamma,
               _isSoc;
    Vecd       _kPath;
    Mat33d     _grid;
    Vector3i   _nGrid;
    MatrixX3i  _gVectors;

  private:
    void _calcKPath();
    const Vector3i get_ngrid( const Mat33d& Acell) const;
    const MatrixX3i gen_gvectors( const Vector3i&  ngrid,      
                                  const int        ikpoint,
                                  const bool       gamma_mode, 
                                  const bool       check_consistency) const;
    const Cubcd gen_ks_wave( const int        ispin,
                             const int        ikpoint,
                             const int        iband,
                                   MatrixX3i  gvectors,
                                   Vector3i   ngrid,
                             const Veccd&     coeff_vec,
                             const int        rescale, 
                             const bool       is_norm) const;
    void save_as_vesta(const Cubcd& phi,
                       const char* header,
                       const char* prefix,
    /* 1=real, 0=im */ const int   mode) const;

}; // end of WAVEHIGH


} // end of namespace

#endif // WAVEHIGH_H
