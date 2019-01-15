#pragma once
#ifndef FFT_H
#define FFT_H


/*
 * FFTW3 Function reference:
 * 
 * Complex DFT
 * 
 *  fftw_plan fftw_plan_dft_1d (int           n,
 *                              fftw_complex* in,   fftw_complex* out,
 *                              int           sign, unsigned      flags);
 * 
 *  fftw_plan fftw_plan_dft_2d (int           n0,   int           n1,
 *                              fftw_complex* in,   fftw_complex* out,
 *                              int           sign, unsigned      flags);
 * 
 *  fftw_plan fftw_plan_dft_3d (int           n0,   int           n1,        int           n2,
 *                              fftw_complex* in,   fftw_complex* out,
 *                              int           sign, unsigned      flags);
 * 
 *  fftw_plan fftw_plan_dft    (int           rank, const int*    n,
 *                              fftw_complex* in,   fftw_complex* out,
 *                              int           sign, unsigned      flags);
 * 
 * 
 * Real DFT
 * 
 *  fftw_plan fftw_plan_dft_r2c_1d(int n,
 *                                 double*  in,   fftw_complex* out,
 *                                 unsigned flags                  );
 * 
 *  fftw_plan fftw_plan_dft_r2c_2d(int      n0,   int           n1,
 *                                 double*  in,   fftw_complex* out,
 *                                 unsigned flags                  );
 * 
 *  fftw_plan fftw_plan_dft_r2c_3d(int      n0,   int           n1, int         n2,
 *                                 double*  in,   fftw_complex* out,
 *                                 unsigned flags                  );
 * 
 *  fftw_plan fftw_plan_dft_r2c   (int      rank, const int*    n,
 *                                 double*  in,   fftw_complex* out,
 *                                 unsigned flags                  );
 * 
 * 
 * Real-to-real DFT
 * 
 *  fftw_plan fftw_plan_r2r_1d(int           n, 
 *                             double*       in,    double*       out,
 *                             fftw_r2r_kind kind,               
 *                             unsigned      flags                  );
 *                                                               
 *  fftw_plan fftw_plan_r2r_2d(int           n0,    int           n1, 
 *                             double*       in,    double*       out,
 *                             fftw_r2r_kind kind0, fftw_r2r_kind kind1, 
 *                             unsigned      flags);
 * 
 *  fftw_plan fftw_plan_r2r_3d(int           n0,     int       n1, int       n2,
 *                             double*       in,     double*   out,
 *                             fftw_r2r_kind kind0,
 *                             fftw_r2r_kind kind1,
 *                             fftw_r2r_kind kind2,
 *                             unsigned      flags);
 * 
 *  fftw_plan fftw_plan_r2r   (int           rank, const int*      n, 
 *                             double*       in,   double*         out, 
 *                             const fftw_r2r_kind*  kind, 
 *                             unsigned      flags);
 * 
 * 
 * For more detailed info: https://www.cnblogs.com/aiguona/p/9407425.html
 */


#include <base.hpp>

namespace ionizing {
  /*
   * using Vecd   = VectorXd;
   * using Veccd  = VectorXcd;
   * using Matd   = Matrix<              double, Dynamic, Dynamic>;
   * using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic>;
   * using Cubd   = Tensor<              double,                3, RowMajor>;
   * using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
   * using Mat33d = Matrix<double,                     3,       3>;
   * using MatX3d = Matrix<double,               Dynamic,       3>;
   */

  Veccd  fft_1d(const Veccd& vec);
  Veccd  fft_1d(const  Vecd& vec);
  Veccd ifft_1d(const Veccd& vec);
  Veccd ifft_1d(const  Vecd& vec);

  Matcd  fft_2d(const Matcd& mat);
  Matcd  fft_2d(const  Matd& mat);
  Matcd ifft_2d(const Matcd& mat);
  Matcd ifft_2d(const  Matd& mat);

  Cubcd  fft_3d(const Cubcd& cub);
  Cubcd  fft_3d(const  Cubd& cub);
  Cubcd ifft_3d(const Cubcd& cub);
  Cubcd ifft_3d(const  Cubd& cub);

  Veccd  rfft_1d(const  Vecd& vec);
  Vecd  irfft_1d(       Veccd vec);
           
  Matcd  rfft_2d(const  Matd& mat);
  Matd  irfft_2d(       Matcd mat);
           
  Cubcd  rfft_3d(const  Cubd& cub);
  Cubd  irfft_3d(       Cubcd cub);
}

#endif // FFT_H
