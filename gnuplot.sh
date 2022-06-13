# p 'SIR.table' u ($1):2 t 'S',\
  '' u ($1):3 t 'E',\
  '' u ($1):4 t 'I',\
  '' u 1:5 t 'R'

# p 'SIR_102400_6.table' every 500 u ($1):2 t 'S',\
  '' every 500 u ($1):3 t 'E',\
  '' every 500 u ($1):4 t 'I',\
  '' every 500 u 1:5 t 'R'

# p 'SIR_102400_4.table' u ($1):2 t 'S',\
  '' u ($1):4 t 'I',\
  '' u 1:5 t 'R'

# p 'SIR_102400_0.2_0.35_4.table' u ($1):($4*102400) t 'I'

# set yr[0.9:1]
# p 'SIR_102400_3_0.35_4.table' u ($1):2 w l lw 3 t 'S',\
  ''                           u ($1):4 w l lw 3 t 'I',\
  ''                           u 1:5    w l lw 3 t 'R'


# p 'SIR_102400_3_0.35_1.table' u ($1):2 w l lw 3 t 'S',\
  ''                           u ($1):4 w l lw 3 t 'I',\
  ''                           u 1:5    w l lw 3 t 'R'


p 'SIR_b4_l6_409600_0.1_0.35_1.table' every 10 u ($1):2 w l lw 3 t 'S',\
  ''                              every 10 u ($1):4 w l lw 3 t 'I',\
  ''                              every 10 u 1:5    w l lw 3 t 'R'
