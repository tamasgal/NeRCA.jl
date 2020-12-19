"""
# KM3NeT Data Definitions v1.3.1-8-gae7dfb2
https://git.km3net.de/common/km3net-dataformat
"""

module Reconstruction
  const JPP_RECONSTRUCTION_TYPE = 4000
  const JMUONBEGIN = 0
  const JMUONPREFIT = 1
  const JMUONSIMPLEX = 2
  const JMUONGANDALF = 3
  const JMUONENERGY = 4
  const JMUONSTART = 5
  const JLINEFIT = 6
  const JMUONEND = 99
  const JSHOWERBEGIN = 100
  const JSHOWERPREFIT = 101
  const JSHOWERPOSITIONFIT = 102
  const JSHOWERCOMPLETEFIT = 103
  const JSHOWER_BJORKEN_Y = 104
  const JSHOWERENERGYPREFIT = 105
  const JSHOWERPOINTSIMPLEX = 106
  const JSHOWERDIRECTIONPREFIT = 107
  const JSHOWEREND = 199
  const DUSJ_RECONSTRUCTION_TYPE = 200
  const DUSJSHOWERBEGIN = 200
  const DUSJSHOWERPREFIT = 201
  const DUSJSHOWERPOSITIONFIT = 202
  const DUSJSHOWERCOMPLETEFIT = 203
  const DUSJSHOWEREND = 299
  const AANET_RECONSTRUCTION_TYPE = 101
  const AASHOWERBEGIN = 300
  const AASHOWERFITPREFIT = 302
  const AASHOWERFITPOSITIONFIT = 303
  const AASHOWERFITDIRECTIONENERGYFIT = 304
  const AASHOWEREND = 399
  const JUSERBEGIN = 1000
  const JMUONVETO = 1001
  const JMUONPATH = 1003
  const JMCEVT = 1004
  const JUSEREND = 1099
  const RECTYPE_UNKNOWN = -1
  const RECSTAGE_UNKNOWN = -1
end