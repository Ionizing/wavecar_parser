#include <wavehigh.hpp>

// #define DEBUG

namespace ionizing {

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  WAVEHIGH::WAVEHIGH(const char* wavecar) : WAVECAR(wavecar) {
#ifdef DEBUG
    std::cout << "WAVEHIGH() created" << std::endl;
#endif
    _calcKPath();
    _nGrid = get_ngrid(getLatticeVectors());

    _gammaHalf = 'z';
    _isGamma = false;
    _isSoc = false;
  } 

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  const Vecd& WAVEHIGH::getKPath() const {
    if (_kPath.size() == 0) {
      std::cerr << " ERROR: _kPath empty" << std::endl;
    }
    
    return _kPath;
  } // end of getKPath()

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  void WAVEHIGH::_calcKPath() {
    const auto& kvec = getKVectors();
    Matrix<double, Dynamic, 3> tmp;
    tmp.resize(kvec.rows() - 1, 3);
    for (int i=0; i!=tmp.rows(); ++i) {
      tmp.row(i) = kvec.row(i + 1).array() - kvec.row(i).array();
    }

    auto tmp2 = tmp * getReciprocalVectors().transpose();
    // std::cout << "\n\n\n" << __LINE__ << " tmp2 = \n" << tmp2 << std::endl;
    _kPath.resize(kvec.rows());
    for (int i=0; i!=kvec.rows(); ++i) {
      if (0 == i) {
        _kPath(i) = 0;
      } else {
        // std::cout << __LINE__ << "tmp2.row(" << i - 1 << ").norm() = " << tmp2.row(i - 1).norm() << std::endl;
        _kPath(i) = _kPath(i - 1) + tmp2.row(i - 1).norm();
      }
    }

    // std::cout << "\n\n\n kpath = \n" << _kPath << std::endl;

  } // end of _calcKPath()

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  void WAVEHIGH::plotBand() const {
    std::ofstream ofs("plotBand.txt");

    if (ofs.fail()) {
      std::cerr << "plotBand.txt open failed!" << std::endl;
      std::abort();
    }

    const auto& band = getBands()[0].real();
    const VectorXd& kpath = getKPath();
    // std::cout << "kpath == \n" << kpath << std::endl;
    // std::cout << "band == \n" << band << std::endl;
    
    MatrixXd res(kpath.rows(), kpath.cols() + band.cols());
    res << kpath, band;
    ofs << res;
    ofs.close();
  } // end of plotBand()

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  const Vector3i WAVEHIGH::get_ngrid(const Mat33d& Acell) const {
    const double encut = getHeader()._enCut;
    Vector3i out;
    for (long i=0; i!=Acell.rows(); ++i) {
      double vec_len = Acell.row(i).norm();
      out(i) = 2 * std::ceil(
            std::sqrt(encut / RY_TO_EV) / (PIx2 / (vec_len / AU_TO_A))
          ) + 1;
    }

#ifdef DEBUG
      std::cout << __FILE__ << __LINE__ << __FUNCTION__ << "ngrid = \n"
        << out << "\n\n\n" << std::endl;
#endif // DEBUG  test passed

    return out;
  } // end of get_ngrid


  const Vector3i& WAVEHIGH::getNGrid() const {
    return _nGrid;
  } // end of getNGrid()

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  const MatrixX3i WAVEHIGH::gen_gvectors(
      const Vector3i& ngrid, const int ikpoint,
      bool gamma_mode, bool check_consistency) const {

    if (ikpoint >= getHeader()._nKpoints or ikpoint < 0) {
      std::cerr << "******** ERROR ********\n" <<
        __FUNCTION__ << ": invalid ikpoint: " << ikpoint << std::endl;
      std::abort();
    }

    ColVecT<Vecd> freqs(3);
    Vecd& freqX = freqs(0);     freqX.resize(ngrid(0));
    Vecd& freqY = freqs(1);     freqY.resize(ngrid(1));
    Vecd& freqZ = freqs(2);     freqZ.resize(ngrid(2));

    for (size_t dim=0; dim!=3; ++dim) { // dim -> dimension: x, y, z
      for (int i=0; i!=ngrid(dim); ++i) {
        // generates: 1, 2, 3, ... ngrid[dim], ..., -3, -2, -1.
        freqs(dim)(i) = (i < ngrid(dim) / 2 + 1) ?  i : i - ngrid(dim);
      }
    }

    // std::cout << __LINE__ << " freqX: " << freqX << std::endl;
    // std::cout << __LINE__ << " freqY: " << freqY << std::endl;
    // std::cout << __LINE__ << " freqZ: " << freqZ << std::endl;

    MatX3d kgrid;
    std::vector<Vector3d> tmpkgrid;
    bool gamma_only = gamma_mode ? true : _isGamma;
    if (gamma_only) {
      if ('z' == _gammaHalf) {
        for (int z=0; z!=ngrid(2); ++z) {
          for (int y=0; y!=ngrid(1); ++y) {
            for (int x=0; x!=ngrid(0); ++x) {
              if ( freqZ(z) > 0 or 
                  (0 == freqZ(z) and freqY(y) > 0) or
                  (0 == freqZ(z) and 0 == freqY(y) and freqX(x) >= 0)) {
                tmpkgrid.push_back({freqX(x), freqY(y), freqZ(z)});
              }
            } // x
          } // y
        } // z
      } else { // 'z' != _gammaHalf
        for (int z=0; z!=ngrid(2); ++z) {
          for (int y=0; y!=ngrid(1); ++y) {
            for (int x=0; x!=ngrid(0); ++x) {
              if ( freqX(x) > 0 or 
                  (0 == freqX(x) and freqY(y) > 0) or
                  (0 == freqX(x) and 0 == freqY(y) and freqZ(z) >= 0)) {
                tmpkgrid.push_back({freqX(x), freqY(y), freqZ(z)});
              }
            } // x
          } // y
        } // z
      } // end if ('z' == _gammaHalf)
    } else { // false == gamma_only
      for (int z=0; z!=ngrid(2); ++z) {
        for (int y=0; y!=ngrid(1); ++y) {
          for (int x=0; x!=ngrid(0); ++x) {
              tmpkgrid.push_back({freqX(x), freqY(y), freqZ(z)});
              // std::cout << __LINE__ << " tmpkgrid[last] = " 
                // << tmpkgrid[tmpkgrid.size() - 1].transpose() << std::endl;
          } // x
        } // y
      } // z
    }

    // constexpr double HBAR2D2ME      = HBAR * HBAR / (2 * M_ELECT);
    const double conversion_constant = HBAR2D2ME;

    MatrixX3i out;
    std::vector<Vector3d> tmpout;

    // std::cout << __LINE__ << ": Bcell = \n" << getReciprocalVectors() << std::endl;

    Vector3d kvec = getKVectors().row(ikpoint);
    for (size_t i=0; i!=tmpkgrid.size(); ++i) {
      double tmp_energy = PIx2 * (
             getReciprocalVectors() * (tmpkgrid[i] + kvec)
          ).norm();
      tmp_energy *= tmp_energy * conversion_constant;
      if (tmp_energy < getHeader()._enCut) {
        // std::cout << __LINE__ << " tmp_energy = " << tmp_energy << std::endl;
        tmpout.push_back(tmpkgrid[i]);
      }
    } // end for i
    
    out.resize(tmpout.size(), 3);
    for (size_t i=0; i!=tmpout.size(); ++i) {
      out.row(i) = tmpout[i].cast<int>();
    }

    if (check_consistency) {
      if (_isSoc) {
        if (out.rows() != getNPlaneWaves()(ikpoint) / 2) {
          fprintf(stderr, "Number of planewaves not consistent with SOC WAVECAR!!\n\
              \tikpoint: %5d, gvector_size: %5ld, _nPlaneWaves(ikpoint): %5d, prod of ngrid: %5d\n", 
            ikpoint, out.rows(), getNPlaneWaves()(ikpoint), ngrid(0) * ngrid(1) * ngrid(2));
          std::abort();
        } 
      } else { // _isSoc
        if (out.rows() != getNPlaneWaves()(ikpoint)) {
          fprintf(stderr, "Number of planewaves not consistent with WAVECAR!!\n\
              \tikpoint: %5d, gvector_size: %5ld, _nPlaneWaves(ikpoint): %5d, prod of ngrid: %5d\n", 
            ikpoint, out.rows(), getNPlaneWaves()(ikpoint), ngrid(0) * ngrid(1) * ngrid(2));
          std::abort();
        }
      } // end if isSoc
    } // end if check_consistency

#ifdef DEBUG
  std::cout << __FUNCTION__ << ": out = \n" << out << std::endl;
#endif
    
    return out;
  } // end of gen_gvectors()

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  const MatrixX3i WAVEHIGH::getGVectors(const int ikpoint) const {
      if (ikpoint >= getHeader()._nKpoints) {
        std::cerr << "******** Invalid ikpoint index, please check. ********" << std::endl;
        std::abort();
      }
    return gen_gvectors(_nGrid, ikpoint, false, true);
  }

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

    const Cubcd WAVEHIGH::gen_ks_wave( const int        ispin,
                                       const int        ikpoint,
                                       const int        iband,
                                             MatrixX3i  gvectors,
                                             Vector3i   ngrid,
                                       const Veccd&     coeff_vec,
                                       const int        rescale, 
                                       const bool       is_norm) const {
    if (!checkIndex(ispin, ikpoint, iband)) {
      std::cerr << "******** CheckIndex failed in " 
        << __FUNCTION__  <<" . Aborting... ********" << std::endl;
      std::abort();
    }

    if (ngrid.size() == 0) {
      ngrid = _nGrid.array() * 2;
    } else {
      // Array<int, 3, 1> tmpsubtract = ngrid - _nGrid;
      const auto difference = (ngrid - _nGrid);
      bool flag = true;
      for (int i=0; i!=ngrid.size(); ++i) {
        if (difference(i) < 0) {
          flag = false;
        }
      }
// #define DEBUG
#ifdef DEBUG
      std::cout << __FUNCTION__ << ": difference = :\n" << difference << std::endl;
#endif //DEBUG
      if (false == flag) {
        std::cerr << "Minimal FFT grid size: " << _nGrid.transpose() << std::endl
          << " while the input ngrid is: " << ngrid.transpose() << std::endl;
        std::abort();
      }
    }

    if (gvectors.size() == 0) {       // TODO: try to prevent copying gvectors.
      gvectors = gen_gvectors(ngrid, ikpoint, false, true);
                  //  ngrid  ikpoint  gamma_mode  check_consistency
    }    

    // gvectors %= 3;
    for (int i=0; i!=3; ++i) {
      int gvec_size = gvectors.rows();
      for (int j=0; j!=gvec_size; ++j) {
        gvectors(j, i) %= ngrid(i);
      }
    }

    double norm_factor = 
        (0 == rescale) ? 
          std::sqrt(ngrid.prod()) : rescale;

    Cubcd phi_ks;  

    if (_isGamma) {
      if (_gammaHalf == 'z') {
        phi_ks.resize(ngrid(0), ngrid(1), ngrid(2) / 2 + 1);
      } else {
        phi_ks.resize(ngrid(0) / 2 + 1, ngrid(1), ngrid(2));
      }
    } else {
      phi_ks.resize(ngrid(0), ngrid(1), ngrid(2) );
    } 

    // std::cerr << __LINE__ << std::endl; std::abort();

    phi_ks.setZero();
    // Cubcd wave_spinor;
    
    for (int i=0; i!=gvectors.rows(); ++i) {
      // gvectors(i, 0) %= ngrid(0);
      // gvectors(i, 1) %= ngrid(1);
      // gvectors(i, 2) %= ngrid(2);
      gvectors(i, 0) = (gvectors(i, 0) < 0) ? ngrid(0) + gvectors(i, 0) : gvectors(i, 0);
      gvectors(i, 1) = (gvectors(i, 1) < 0) ? ngrid(1) + gvectors(i, 1) : gvectors(i, 1);
      gvectors(i, 2) = (gvectors(i, 2) < 0) ? ngrid(2) + gvectors(i, 2) : gvectors(i, 2);
    }

    // std::cout << __FUNCTION__ << ": gvectors = \n" << gvectors << std::endl;

    if (_isSoc) {
      std::cerr << "SOC WAVECAR gen_ks_wave not supported yet, impl later..." << std::endl;
      std::abort();

/*
 *       Veccd tmp_coeff;
 *       if (coeff_vec.size() != 0) {
 *         tmp_coeff = coeff_vec;
 *       } else {
 *         tmp_coeff = getBandCoeff(ispin, ikpoint, iband, is_norm);
 *       }
 *       
 *       const int n_plane_waves = tmp_coeff.size() / 2;
 * 
 *       
 * 
 *       // spinor up
 *       for (int i=0; i!=n_plane_waves; ++i) {
 *         const int idx = gvectors(i, 0);
 *         const int idy = gvectors(i, 1);
 *         const int idz = gvectors(i, 2);
 *         phi_ks(idx, idy, idz) = tmp_coeff(i);
 *       }
 * 
 *       for (int i=n_plane_waves; i!=tmp_coeff.size(); ++i) {
 *         const int idx = gvectors(i, 0);
 *         const int idy = gvectors(i, 1);
 *         const int idz = gvectors(i, 2);
 *         phi_ks(idx, idy, idz) = tmp_coeff(i);
 *       }
 */

    } else { // not SOC


    // std::cerr << __LINE__ << std::endl; std::abort();

      Veccd tmp_coeff;
      if (coeff_vec.size() != 0) {
        tmp_coeff = coeff_vec;
      } else {
        tmp_coeff = getBandCoeff(ispin, ikpoint, iband, is_norm);
      }

#ifdef DEBUG
      { // debug
        std::ofstream ofs("tmp_coeff.txt");
        ofs << tmp_coeff << std::endl;
        ofs.close();
      }
#endif
      
      const int n_plane_waves = tmp_coeff.size();
      for (int i=0; i!=n_plane_waves; ++i) {
        const int idx = gvectors(i, 0);
        const int idy = gvectors(i, 1);
        const int idz = gvectors(i, 2);
        phi_ks(idx, idy, idz) = tmp_coeff(i);
      }

#ifdef DEBUG
        { // debug
          std::ofstream ofs("phi_ks_unnormed.txt");
          ofs << phi_ks.chip(0, 0) << std::endl;
          ofs.close();
        }
#endif

      if (_isGamma) {
        std::cerr << "gamma mode not finished yet, impl later." << std::endl;
        std::abort();
        
/*
 *         if ('z' == _gammaHalf) {
 *           for (int i=0; i!=ngrid(0); ++i) {
 *             for (int j=ngrid(1)/2 + 1; j!=ngrid(1); ++j) {
 *               const int freqX = (i < ngrid(0)/2 + 1) ? i : i - ngrid(0);
 *               const int freqY = j - ngrid(1);
 *               if ((freqY > 0) or (freqY == 0 and freqX >= 0)) {
 *                 continue;
 *               } else {
 *                 const int idx = i % ngrid(0);
 *                 const int idy = j % ngrid(1);
 *                 phi_ks(idx, idy, 0) = 
 *                     std::conj(phi_ks(ngrid(0) - idx, ngrid(1) - idy, 0));
 *               }
 *             } // end for j
 *           } // end for i
 *           
 *           phi_ks /= std::sqrt(2.0);
 *           phi_ks(0, 0, 0) *= std::sqrt(2.0);
 *           return ifft_3d(  )
 * 
 * 
 *         } else {  // gammaHalf != 'z'
 * 
 * 
 * 
 *         }
 */
      } else {
        // Most probable condition
#ifdef DEBUG
        { // DEBUG
          std::cout << __FILE__ << __LINE__ << "  norm_factor = " << norm_factor << std::endl;
        }
#endif


        // use sequent access to accelerate
        for (int i=0; i!=phi_ks.size(); ++i) {
          *(phi_ks.data() + i) *= norm_factor;
        }
        // slower version
        /*
         * for (int i=0; i!=phi_ks.dimension(0); ++i) {
         *   for (int j=0; j!=phi_ks.dimension(1); ++j) {
         *     for (int k=0; k!=phi_ks.dimension(2); ++k) {
         *       phi_ks(i, j, k) *= norm_factor;
         *     }
         *   }
         * }
         */

#ifdef DEBUG
        { // debug
          std::ofstream ofs("phi_ks_normed.txt");
          ofs << phi_ks.chip(0, 0) << std::endl;
          ofs.close();
        }
#endif

        auto out = ifft_3d(phi_ks);

#ifdef DEBUG
        { // debug
          std::ofstream ofs("phi_ks_ifft.txt");
          ofs << out.chip(0, 0) << std::endl;
          ofs.close();
        }
#endif

        return out;

      }
      

    }  // Not SOC
    
    std::cerr << "Unhandled condition in " << __FUNCTION__ << std::endl;
    std::abort(); 
  }

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

    const Cubcd WAVEHIGH::getKSWave(const int       ispin,
                                    const int       ikpoint,
                                    const int       iband,
                                          MatrixX3i gvectors,
                                          Vector3i  ngrid,
                                    const Veccd&    coeff_vec,
                                    const int       rescale,
                                    const bool      is_norm) const {
      if (ispin >= getInfo()._nSpin or
          ikpoint >= getHeader()._nKpoints or
          iband >= getHeader()._nBands) {
        std::cerr << "******** Invalid ispin, ikpoint, iband index, please check. ********" << std::endl;
        std::abort();
      }
      return gen_ks_wave(ispin, ikpoint, iband, gvectors,
          ngrid, coeff_vec, rescale, is_norm);

    }


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

    void WAVEHIGH::saveAsVesta(const Cubcd& phi,
                               const char*  POSCAR,
                               const char*  prefix,
                               const bool   is_real) const {

      char file_name[256];
      if (strcmp(POSCAR, "") == 0) {
        strcpy(file_name, "POSCAR");
      } else {
        strcpy(file_name, POSCAR);
      }

      if (strlen(file_name) >= 255) {
        std::cerr << "POSCAR file name too long." << std::endl;
        std::abort();
      }
      std::ifstream poscar(file_name);
      if (poscar.fail()) {
        std::cerr << "POSCAR: " << file_name <<" open failed !" << std::endl;
        std::abort();
      }

      /***** Generating Header *****/
      std::string header;
      while (!poscar.eof()) {
        std::string line;       
        std::getline(poscar, line);
        if (trim_copy(line).size() == 0) {
          header += "\n";
          break;
        }
        line += "\n";
        header += line;
      }
      poscar.close();

      // write real part
      save_as_vesta(phi, header.c_str(), prefix, 1);

#ifdef DEBUG
      { // DEBUG
        std::cout << __FILE__ << __LINE__ 
          << " _isGamma = " << std::boolalpha << _isGamma
          << " is_real = " << is_real << std::endl;
      }
#endif

      if (!(_isGamma or is_real)) {
        // write image part
        save_as_vesta(phi, header.c_str(), prefix, 0);
      }

      return;
    }

  

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

    void WAVEHIGH::save_as_vesta(const Cubcd& phi,
                                 const char* header,
                                 const char* prefix,
   /* 1 = real, 0 = imagine  */  const int mode) const {
      char out_file_name[256];
      if (strlen(prefix) >= 255-7) {
        std::cerr << "Prefix too long." << std::endl;
        std::abort();
      }

      if (0 == mode) {
        sprintf(out_file_name, "%s%s", prefix, "_i.vasp");
      } else {
        sprintf(out_file_name, "%s%s", prefix, "_r.vasp");
      }

      std::ofstream out(out_file_name);
      out << header;


      const int nx = phi.dimension(0),
                ny = phi.dimension(1),
                nz = phi.dimension(2);
      char line[1<<20]; // buffer size == 10M
      sprintf(line, "%5d%5d%5d\n", nx, ny, nz);
      out << line;
      strcpy(line, "");

      if (0 == mode) {
        int nwrite = 0;
        for (int k=0; k!=nz; ++k) {
          for (int j=0; j!=ny; ++j) {
            for (int i=0; i!=nx; ++i) {
              nwrite += 1;
              sprintf(line + strlen(line), "%16.8E ", phi(i, j, k).imag());
              if (nwrite % 10 == 0) {
                strcat(line, "\n");
                out << line;
                strcpy(line, "");
              }
            } // i
          } // j
        } // k
      } else {
        int nwrite = 0;
        for (int k=0; k!=nz; ++k) {
          for (int j=0; j!=ny; ++j) {
            for (int i=0; i!=nx; ++i) {
              nwrite += 1;
              sprintf(line + strlen(line), "%16.8E ", phi(i, j, k).real());
              if (nwrite % 10 == 0) {
                strcat(line, "\n");
                out << line;
                strcpy(line, "");
              }
            } // i
          } // j
        } // k
      }

      out.close();
      
      return;
    }


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

  void WAVEHIGH::plotWave(const int ispin,
                          const int ikpoint,
                          const int iband,
                          const char* poscar,
                          const char* prefix,
                          const bool is_real) const {
    
    if (ispin >= getInfo()._nSpin or
        ikpoint >= getHeader()._nKpoints or
        iband >= getHeader()._nBands) {
      std::cerr << "******** Invalid ispin, ikpoint, iband index, please check. ********" << std::endl;
      std::abort();
    }

    const Vector3i ngrid = getNGrid().array() * 2;
    const MatrixX3i gvecs = getGVectors(ikpoint);
    const Veccd& coeff_vec = getBandCoeff(ispin, ikpoint, iband, true);
    const auto phi = getKSWave(ispin, ikpoint, iband, 
        gvecs, ngrid, coeff_vec, 0, true);

    using std::cout;
    using std::setw;
    using std::setfill;
    using std::endl;
    std::stringstream ss;
    ss << prefix << "_" << setfill('0') << setw(2) << ispin << "_" 
       << setw(2) << ikpoint << "_" << setw(2) << iband;
    saveAsVesta(phi, poscar, ss.str().c_str(), is_real);

    cout << "plotWave at spin " << setw(2) << ispin
      << " kpoint " << setw(2) << ikpoint
      << " band " << setw(2) << iband << " finished." << endl;
  }




/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

} // end of namesace ionizing
