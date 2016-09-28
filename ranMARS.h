/*************************************************

      ranMARS:
      MARSAGLIA pseudo random number generator

*************************************************/
#ifndef _RANMARS_H_
#define _RANMARS_H_

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
  
#define MARS_FIELD_SIZE 98

  
  void init_RAN( int nA1, int nA2, int nA3, int nB1 );
  void save_RAN( char *fileName );
  void restore_RAN( char *fileName ); 
  double RAN01( void );

#ifdef __cplusplus
}
#endif

#endif // _RANMARS_H_


