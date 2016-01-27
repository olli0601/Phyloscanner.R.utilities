#include <cstdarg>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <cmath>
#include <Rmath.h>
#include "hivc.h"

// ERROR HANDLING
static inline void postIfError(const char* format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
	fflush(stdout);
 	std::exit(EXIT_FAILURE);
 	//throw std::runtime_error(nabcGlobals::BUFFER);
}
#define FAIL_ON(condition, message) if (condition){ postIfError(message, condition); }/**< macro to throw error messages */

// ASSUMED FROM APE
//	tab_trans[65] = 0x88; /* A */
//	tab_trans[71] = 0x48; /* G */
//	tab_trans[67] = 0x28; /* C */
// 	tab_trans[84] = 0x18; /* T */
// 	tab_trans[82] = 0xc0; /* R */
// 	tab_trans[77] = 0xa0; /* M */
// 	tab_trans[87] = 0x90; /* W */
// 	tab_trans[83] = 0x60; /* S */
// 	tab_trans[75] = 0x50; /* K */
// 	tab_trans[89] = 0x30; /* Y */
// 	tab_trans[86] = 0xe0; /* V */
// 	tab_trans[72] = 0xb0; /* H */
// 	tab_trans[68] = 0xd0; /* D */
//  	tab_trans[66] = 0x70; /* B */
// 	tab_trans[78] = 0xf0; /* N */
//
//	tab_trans[97] = 0x88; /* a */
//	tab_trans[103] = 0x48; /* g */
//	tab_trans[99] = 0x28; /* c */
// 	tab_trans[116] = 0x18; /* t */
// 	tab_trans[114] = 0xc0; /* r */
// 	tab_trans[109] = 0xa0; /* m */
// 	tab_trans[119] = 0x90; /* w */
// 	tab_trans[115] = 0x60; /* s */
// 	tab_trans[107] = 0x50; /* k */
// 	tab_trans[121] = 0x30; /* y */
// 	tab_trans[118] = 0xe0; /* v */
// 	tab_trans[104] = 0xb0; /* h */
// 	tab_trans[100] = 0xd0; /* d */
//  	tab_trans[98] = 0x70; /* b */
// 	tab_trans[110] = 0xf0; /* n */
//
//  	tab_trans[45] = 0x04; /* - */
//  	tab_trans[63] = 0x02; /* ? */

/*
------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
A            A           Adenine
C            C           Cytosine
G            G           Guanine
T            T           Thymine
U            U           Uracil
M          A or C
R          A or G
W          A or T
S          C or G
Y          C or T
K          G or T
V        A or C or G
H        A or C or T
D        A or G or T
B        C or G or T
X      G or A or T or C
N      G or A or T or C
*/


/* these return 1 if the base is Base, 0 otherwise */
#define IsA(a) (a == 0x88)
#define IsG(a) (a == 0x48)
#define IsC(a) (a == 0x28)
#define IsT(a) (a == 0x18)
#define IsR(a) (a == 0xc0)
#define IsM(a) (a == 0xa0)
#define IsW(a) (a == 0x90)
#define IsS(a) (a == 0x60)
#define IsK(a) (a == 0x50)
#define IsY(a) (a == 0x30)
#define IsV(a) (a == 0xe0)
#define IsH(a) (a == 0xb0)
#define IsD(a) (a == 0xd0)
#define IsB(a) (a == 0x70)
#define IsN(a) (a == 0xf0)
#define IsGap(a) (a == 0x04)
#define IsNA(a) (a == 0x02)

#define MatchesA(b)	(IsA(b) || IsM(b) || IsR(b) || IsW(b) || IsV(b) || IsH(b) || IsD(b) || IsN(b))
#define MatchesC(b)	(IsC(b) || IsM(b) || IsS(b) || IsY(b) || IsV(b) || IsH(b) || IsB(b) || IsN(b))
#define MatchesG(b)	(IsG(b) || IsR(b) || IsS(b) || IsK(b) || IsV(b) || IsD(b) || IsB(b) || IsN(b))
#define MatchesT(b)	(IsT(b) || IsW(b) || IsY(b) || IsK(b) || IsH(b) || IsD(b) || IsB(b) || IsN(b))

#define IsPairA(a,b)	(MatchesA(a) && MatchesA(b))
#define IsPairC(a,b)	(MatchesC(a) && MatchesC(b))
#define IsPairG(a,b)	(MatchesG(a) && MatchesG(b))
#define IsPairT(a,b)	(MatchesT(a) && MatchesT(b))

#define IsPair(a,b)		(IsPairA(a,b) || IsPairC(a,b) || IsPairG(a,b) || IsPairT(a,b))

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* these return 1 if the base is known or ambiguous, 0 otherwise */
#define IsKnownOrAmbiguousBase(a) (!IsGap(a) & !IsNA(a))


void hivc_printdna(unsigned char *x, int *n)
{
	int i=0;
	for(; i<*n;i++)
	{
		std::cout<<KnownBase(x[i])<<'\t'<<IsA(x[i])<<'\t'<<IsC(x[i])<<'\t'<<IsG(x[i])<<'\t'<<IsT(x[i])<<'\t'<<IsM(x[i])<<'\t'<<IsR(x[i])<<'\t'<<IsW(x[i])<<'\t'<<IsS(x[i])<<'\t'<<IsY(x[i])<<'\t'<<IsK(x[i])<<'\t'<<IsV(x[i])<<'\t'<<IsH(x[i])<<'\t'<<IsD(x[i])<<'\t'<<IsB(x[i])<<'\t'<<IsN(x[i])<< std::endl;
	}
}

void hivc_dist_ambiguous_dna(unsigned char *x1, unsigned char *x2, int *n, double *ans)
{
	int i= *n, m1=0, m2=0;
	unsigned char *y1= x1, *y2= x2;

	for(; 	i--; 	y1++, y2++)
	{
		if(IsPair(*y1,*y2))	m1++;
		//if(IsPair(*y1,*y2))	std::cout<<"Fnd Pair at pos"<<(*n-i+1)<<std::endl;
		//if(MatchesC(*y1))	std::cout<<"Fnd y1 MatchC at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(MatchesC(*y2))	std::cout<<"Fnd y2 MatchC at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(MatchesT(*y2))	std::cout<<"Fnd y2 MatchT at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(IsPairA(*y1,*y2))	std::cout<<"Fnd PairA at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(IsPairC(*y1,*y2))	std::cout<<"Fnd PairC at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(IsPairG(*y1,*y2))	std::cout<<"Fnd PairG at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(IsPairT(*y1,*y2))	std::cout<<"Fnd PairT at pos"<<(*n-i-1)<<" "<<m1<<std::endl;
		//if(IsA(*y1))	std::cout<<"Fnd IsA at pos"<<(*n-i-1)<<std::endl;
		if(IsKnownOrAmbiguousBase(*y1) && IsKnownOrAmbiguousBase(*y2))	m2++;
	}
	//std::cout<<"ans "<<m1<<"/"<<m2<<std::endl;
	*ans= static_cast<double>(m1)/m2;
}

SEXP hivc_clu_mintransmissioncascade(SEXP inbrlv)
{
	FAIL_ON(! Rf_isVector(inbrlv) ,"hivc_clu_mintransmissioncascade: brlv is not vector %c");
	SEXP ans;
	PROTECT(ans= Rf_allocVector(REALSXP, 1));

	int i, j, m, nbrl=Rf_length(inbrlv), nchange=0, nequal=0, verbose=0;
	double tmp= std::sqrt( nbrl*2+0.25 ) + 0.5;
	FAIL_ON( std::floor(tmp)!=tmp, "hivc_clu_mintransmissioncascade: brlv is not vector of upper triangular matrix %c");
	int nc= CAST(int, tmp);
	double *xd= NULL, *yd= NULL, *brlv=NULL;
	double **xbrlv=NULL, **ybrlv= NULL, **triangle=NULL;

	if(nc==1)
		*REAL(ans)= 0;
	else if(nc==2)
		*REAL(ans)= *REAL(inbrlv);
	else
	{
		//make internal copy of SEXP
		brlv= NEW_ARY(double,nbrl);
		for(i= nbrl, xd= REAL(inbrlv), yd=brlv;  	i; 		i--, *yd++= *xd++);
		//set up array for triangular comparison
		triangle= NEW_ARY(double*,3);
		//set up fast indexing
		xbrlv= NEW_ARY(double*,nc-1);
		for(i= nc-1, xd= brlv, ybrlv=xbrlv;  	i; 		i--)
		{
			*ybrlv++= xd;
			xd+= i;
		}
		if(verbose)
		{
			std::cout<<"\nncol is "<<nc<<std::endl;
			std::cout<<"\ninput"<<std::endl;
			for(i=0; 	i<nc;  i++)
			{
				std::cout<<'\n';
				for(j=0;		j<=i; 		j++ )
					std::cout<< 0 <<'\t';
				for(j=0;		j<(nc-i-1); 		j++ )
					std::cout<< *(xbrlv[i]+j) <<'\t';
			}
			std::cout<<std::endl;
			//FAIL_ON(1,"stophere");
		}

		for(m=nc-1;		m--;		)//need at most n-1 rounds of transitive closure
		{
			//if all equal, can abort
			for(i=nbrl, tmp= *(xd= brlv), nequal=0; 	i--;  nequal++)
				if(tmp!=*xd++)
					break;
			if(verbose)
				std::cout<<"\nnequal "<<nequal<<std::endl;
			if(nequal==nbrl)
				break;
			//one round of transitive closure
			for(i=0, nchange=0; 	i<(nc-2);		i++)
				for(j=1;		j<(nc-i-1);		j++)
				{
					/*  triangle coding is   0	1
					 * 							2
					 */
					triangle[0]= xbrlv[i];
					triangle[1]= xbrlv[i]+j;
					triangle[2]= xbrlv[i+1]+j-1;

					if(verbose)
						std::cout<<'\n'<<(**triangle)<<'\t'<<(**(triangle+1))<<'\t'<<(**(triangle+2));
					if( **triangle<=**(triangle+1))
					{
						if(			**(triangle+1)<=**(triangle+2))// 0,1,2		ie 2<- 1
						{
							if(**(triangle+2)!=**(triangle+1)) 	nchange++;
							**(triangle+2)= **(triangle+1);
						}
						else if(	**(triangle+1)>**(triangle+2)  && 	**triangle<**(triangle+2))// 0,2,1	 ie 1<- 2
						{
							if(**(triangle+1)!=**(triangle+2)) 	nchange++;
							**(triangle+1)= **(triangle+2);
						}
						else	// 2,0,1 	ie	1<- 0
						{
							if(**(triangle+1)!=**triangle) 		nchange++;
							**(triangle+1)= **triangle;
						}
					}
					else
					{
						if(			**(triangle+1)>**(triangle+2))// 2,1,0 	ie	0<- 1
						{
							if(**triangle!=**(triangle+1)) 	nchange++;
							**triangle= **(triangle+1);
						}
						else if(	**(triangle+1)<=**(triangle+2)  && 	**triangle<=**(triangle+2))// 1,0,2		ie 2<- 0
						{
							if(**(triangle+2)!=**triangle) 	nchange++;
							**(triangle+2)= **triangle;
						}
						else	// 1,2,0	ie 0<- 2
						{
							if(**triangle!=**(triangle+2)) 	nchange++;
							**triangle= **(triangle+2);
						}
					}
					if(verbose)
						std::cout<<'\n'<<(**triangle)<<'\t'<<(**(triangle+1))<<'\t'<<(**(triangle+2));

				}
			if(verbose)
			{
				std::cout<<"\ntriangularization step "<<m<<std::endl;
				for(i=0; 	i<nc;  i++)
				{
					std::cout<<'\n';
					for(j=0;		j<=i; 		j++ )
						std::cout<< 0 <<'\t';
					for(j=0;		j<(nc-i-1); 		j++ )
						std::cout<< *(xbrlv[i]+j) <<'\t';
				}
				std::cout<<std::endl;
				std::cout<<"\nnchange "<<nchange<<std::endl;
			}
			//if no change, can abort
			if(nchange==0)
				break;
		}
		//compute max brl over upper triangular
		for(i=nbrl, tmp= *(xd= brlv); 	i--;  xd++)
			if(tmp<*xd)
				tmp= *xd;
		*REAL(ans)= tmp;
		DELETE(brlv);
		DELETE(triangle);
		DELETE(xbrlv);
	}
	UNPROTECT(1);
	return ans;
}
