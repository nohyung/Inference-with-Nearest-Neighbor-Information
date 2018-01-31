#include "mex.h"


void getNChooseKdouble(double* nchooseks, int n)
{
    double* curarray;
    int icnt, icnt1, icnt2;
//     printf("n[%d]", n);
    
    curarray = (double*)malloc((n + 1)*sizeof(double));
    curarray[0] = 1;
    nchooseks[0] = curarray[0];     // moved location for n == 0
    for (icnt = 1; icnt < n + 1; icnt++)
        curarray[icnt] = 0;
    for (icnt1 = 0; icnt1 < n; icnt1++)
    {
 //       nchooseks[0] = curarray[0]; // do not work when n == 0
        for (icnt2 = 1; icnt2 < icnt1 + 2; icnt2++)
            nchooseks[icnt2] = curarray[icnt2 - 1] + curarray[icnt2];
        if (icnt1 != n - 1) {
            for (icnt2 = 0; icnt2 < icnt1 + 2; icnt2++)
                curarray[icnt2] = nchooseks[icnt2];
        }
    }
    free(curarray);
}

void mexFunction(int ds, mxArray *d[], int ss, const mxArray *s[])
{
    double *nnum, *Postrs;
    int datanum;


//     printf("ds[%d], ss[%d]\n", ds, ss);
    if (ds != 1 || ss != 1 ) mexErrMsgTxt("usage: nchoosenval = NChooseKdouble(nnum)");

    // input data
    nnum = mxGetPr(s[0]);

    // output data
    datanum = *nnum + 1;
    d[0] =  mxCreateDoubleMatrix(1, datanum, mxREAL);
    Postrs = mxGetPr(d[0]);

    getNChooseKdouble(Postrs, (int)*nnum);

// //     ////////////////////////
// //     icnt = 32;
// // //     icnt = 100;
// // //     nchooseks = (int*)malloc((icnt + 1)*sizeof(int));
// // //     getNChooseK(nchooseks, icnt);
// //     nchooseks = (double*)malloc((icnt + 1)*sizeof(double));
// //     getNChooseKdouble(nchooseks, icnt);
// //     for (icnt1 = 0; icnt1 < icnt + 1; icnt1++)
// // //         printf("nchooseks[%d] ", nchooseks[icnt1]);
// //         printf("nchooseks[%f] ", nchooseks[icnt1]);
// //     printf("\n");
// //     free(nchooseks);
// //     ////////////////////////
// 
// 
//     for (icnt = 0; icnt < datanum; icnt++)
//     {
//         curMinClass = nextNeighborClass[icnt] - 1;
//         d1 = distPerClasses[icnt + datanum*curMinClass];
//         Postrs[icnt] = 1;
//         for (iclass = 0; iclass < classnum; iclass++)
//         {
//             if (iclass == curMinClass)
//                 continue;
//             d2 = distPerClasses[icnt + datanum*iclass];
//             PosteriorOneClasss = 0;
//             nchooseks = (double*)malloc((curKnums[icnt + datanum*iclass] + curKnums[icnt + datanum*curMinClass] + 2)*sizeof(double));
//             getNChooseKdouble(nchooseks, (int)(curKnums[icnt + datanum*iclass] + curKnums[icnt + datanum*curMinClass] + 1));
//             for (isum1 = 0; isum1 < curKnums[icnt + datanum*curMinClass] + 1; isum1++)
//             {
// //                 printf("nchooseks[isum1] %f\n", nchooseks[isum1]);
// //                 printf("%f\n", isum1*(log(1 + pow(d1/ *bval, featurenum)) - log(2 + pow(d1/ *bval, featurenum) + pow(d2/ *bval, featurenum))));
// //                 printf("%f\n", (curKnums[icnt + datanum*iclass] + curKnums[icnt + datanum*curMinClass] + 1 - isum1)*(log(1 + pow(d2/ *bval, featurenum)) - log(2 + pow(d1/ *bval, featurenum) + pow(d2/ *bval, featurenum))));
// //                 printf("%f\n", log(1 + pow(d1/ *bval, featurenum)));
// //                 printf("%f\n", pow(d1/ *bval, featurenum));
// //                 printf("%f %f\n", d1, *bval);
// //                 printf("%f \n", *bval);
// //                 printf("%f\n", log(2 + pow(d1/ *bval, featurenum) + pow(d2/ *bval, featurenum)));
//                 PosteriorOneClasss += nchooseks[isum1]*
//                         exp(isum1*(log(1 + pow(d1/ *bval, featurenum)) - log(2 + pow(d1/ *bval, featurenum) + pow(d2/ *bval, featurenum))) +
//                         (curKnums[icnt + datanum*iclass] + curKnums[icnt + datanum*curMinClass] + 1 - isum1)*(log(1 + pow(d2/ *bval, featurenum)) - log(2 + pow(d1/ *bval, featurenum) + pow(d2/ *bval, featurenum))));
//             }
// 
//             free(nchooseks);
//             Postrs[icnt] *= PosteriorOneClasss;
//         }
//     }
// //     if (ss != 5)  free(bval);
// 
//     
// //     for (icnt = 0; icnt < datanum; icnt++)
// //         printf("[%f] ", Postrs[icnt]);
// //     printf("\n");
}
