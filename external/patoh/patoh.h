/*
  ---------------------------------------------------------------------
  This file is part of PaToH (Partitioning Tools for Hypergraphs)
  Copyright (c) 1996-   Umit V. Catalyurek
  ---------------------------------------------------------------------
  For license info, please see the README file in the main directory.
  ---------------------------------------------------------------------

  -------------------------------------------------------------
  File    : patoh.h
  Author  : Umit V. Catalyurek
  Date    : June 2, 2000
  -------------------------------------------------------------
  Description: PaToH v3.2 Library Interface
  -------------------------------------------------------------
*/
#ifndef _PATOH_V3_LIB_H_
#define _PATOH_V3_LIB_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ======================= CONSTANTS ======================= */

/* ------------------- Cut Metrics -------------------*/
/* Conectivity-1 cut definition*/
#define PATOH_CONPART                   1
/* netcut definition */
#define PATOH_CUTPART                   2

/* ------------------- Parameter Initialization -------------------*/
#define PATOH_SUGPARAM_DEFAULT          0
#define PATOH_SUGPARAM_SPEED            1
#define PATOH_SUGPARAM_QUALITY          2

/* ------------------- Balance  -------------------*/
#define PATOH_BALANCE_STRICT            0
#define PATOH_BALANCE_ADAPTIVE          1
#define PATOH_BALANCE_RELAXED           2

/* ------------------- Coarsening  ------------------- */
/* -------- Visit Order -------------- */
#define PATOH_VO_CONT                   0
#define PATOH_VO_RAND                   1
#define PATOH_VO_SCD                    2
#define PATOH_VO_SMAXNS                 3
#define PATOH_VO_SMINNS                 4
#define PATOH_VO_SMINNSSUM              5
#define PATOH_VO_SWEEP                  6



/* -------- Matching -------------- */
#define PATOH_CRS_FF                    0
#define PATOH_CRS_HCM                   1
#define PATOH_CRS_PHCM                  2
#define PATOH_CRS_MANDIS                3
#define PATOH_CRS_AVGDIS                4
#define PATOH_CRS_CANBERA               5
#define PATOH_CRS_ABS                   6
#define PATOH_CRS_GCM                   7
#define PATOH_CRS_SHCM                  8
/* --------- Agglomeratives ---------*/
#define PATOH_CRS_HCC                   9
#define PATOH_CRS_HPC                   10
#define PATOH_CRS_ABSHCC                11
#define PATOH_CRS_ABSHPC                12
#define PATOH_CRS_CONC                  13
#define PATOH_CRS_GCC                   14
#define PATOH_CRS_SHCC                  15
#define PATOH_CRS_FIRST_NET_MATCH       (PATOH_CRS_SHCC+1)
/* --------------- Net Base Agglomeratives -------------- */
#define PATOH_CRS_NC                    PATOH_CRS_FIRST_NET_MATCH
#define PATOH_CRS_MNC                   (PATOH_CRS_FIRST_NET_MATCH+1)




/* ------------------- Inital Partitionin ------------------- */
#define PATOH_IPA_GHGP                  1
#define PATOH_IPA_AGGMATCH              2
#define PATOH_IPA_BF                    3
#define PATOH_IPA_BINPACK               4
#define PATOH_IPA_RANDOM1               5
#define PATOH_IPA_RANDOM2               6
#define PATOH_IPA_RANDOM3               7
#define PATOH_IPA_GHG_MAXPIN            8
#define PATOH_IPA_GHG_MAXNET            9
#define PATOH_IPA_GHG_MAXPOSGAIN        10
#define PATOH_IPA_COMP_GHGP             11
#define PATOH_IPA_GREEDY_COMP_GHGP      12
#define PATOH_IPA_ALL                   13



/* ------------------- Refinement ------------------- */
#define PATOH_REFALG_NONE               0
#define PATOH_REFALG_T_BFM              1
#define PATOH_REFALG_T_FM               2
#define PATOH_REFALG_D_BFM              3
#define PATOH_REFALG_D_FM               4
#define PATOH_REFALG_BKL                5
#define PATOH_REFALG_KL                 6
#define PATOH_REFALG_MLG_BFM            7
#define PATOH_REFALG_MLG_FM             8
#define PATOH_REFALG_BFMKL              9
#define PATOH_REFALG_FMKL               10


/* ------------------- Output Detail ------------------- */
#define PATOH_OD_ONLYRESTIME            -1
#define PATOH_OD_NONE                   0
#define PATOH_OD_LOW                    1
#define PATOH_OD_MEDIUM                 2
#define PATOH_OD_HIGH                   3



/* ======================= TYPES ======================= */


    typedef struct
    {
        /* ================ Miscellaneous Parameters ================ */
        /* ------- general parameters ------------*/
        int             cuttype;
        int             _k;
        int             outputdetail;        /* PATOH_OD_... */
        int             seed;
        int             doinitperm;

        /* --------- net discard  parameters -----------*/
        int             bisec_fixednetsizetrsh;
        float           bisec_netsizetrsh;
        int             bisec_partmultnetsizetrsh;

        /*----------------- V-cycle parameter ----------------*/
        int             bigVcycle,
                        smallVcycle,
                        usesamematchinginVcycles;

        /* ---------- use heap/bucket parameters ---------- */
        int             usebucket;
        int             maxcellinheap;
        int             heapchk_mul;
        int             heapchk_div;

        /* ----------------- Memory Allocation Parameters ------------------*/
        int             MemMul_CellNet,
                        MemMul_Pins,
                        MemMul_General;


        /* ================ Coarsening Parameters ================ */
        int             crs_VisitOrder;      /* PATOH_VO_... */
        int             crs_alg;             /* PATOH_CRS_... */
        int             crs_coarsento,
                        crs_coarsentokmult,
                        crs_coarsenper;
        float           crs_maxallowedcellwmult;
        int             crs_idenafter;
        int             crs_iden_netsizetrh;
        int             crs_useafter, crs_useafteralg;

        /* ================ Initial Partitioning Parameters ================ */
        /*--- both init part & refinement -----*/
        int             nofinstances;

        /* -------------- initial partitioning parameters ---------------*/
        int             initp_alg;           /* PATOH_IPA_... */
        int             initp_runno;
        int             initp_ghg_trybalance;
        int             initp_refalg;        /* PATOH_REFALG_... */

        /* ================ Refinement Parameters ================ */
        int             ref_alg;             /* PATOH_REFALG_... */
        int             ref_useafter, ref_useafteralg;
        int             ref_passcnt,
                        ref_maxnegmove;
        float           ref_maxnegmovemult;
        int             ref_dynamiclockcnt;
        float           ref_slow_uncoarsening;

        /* -------------- imbalance parameters--------------*/
        int             balance;
        double          init_imbal,
                        final_imbal,
                        fast_initbal_mult;
        float           init_sol_discard_mult,
                        final_sol_discard_mult;

/* some utility parameters that might be useful in user's application; not used by PaToH (other than parsing them in Process_Arguments) */
        char            allargs[8192];
        char            inputfilename[512];
        int             noofrun;
        int             writepartinfo;


    } PaToH_Parameters;

    typedef PaToH_Parameters             *PPaToH_Parameters;


/* ======================= FUNCTION PROTOTYPES ======================= */

/* --- inititalization and memory functions --- */
    int PaToH_Initialize_Parameters(PPaToH_Parameters pargs, int cuttype,
                                    int SuggestByProblemType);


    /* Function should return 0 if it does NOT processes the argument, any other integer means it DID process the argument */
    typedef int (*SingleArgumentCallBack)(PPaToH_Parameters pargs, char *variable, char *value, int ivalue, float fvalue);

    int PaToH_Process_Arguments(PPaToH_Parameters pargs, int startargc, int argc, char *argv[], SingleArgumentCallBack func);



    int PaToH_Alloc(PPaToH_Parameters pargs, int _c, int _n, int _nconst,
                    int *cwghts, int *nwghts, int *xpins, int *pins);

    int PaToH_Free(void);

/* --- partition --- */

    /* Input: pargs, _c, _n, cwghts, nwghts, xpins, pins
       Ouput: partvec, partweights and cut */
    /* Deprecated: use PaToH_Part instead */
    int PaToH_Partition(PPaToH_Parameters pargs, int _c, int _n, int *cwghts,
                        int *nwghts, int *xpins, int *pins, int *partvec,
                        int *partweights, int *cut);

    /* Input: pargs, _c, _n, cwghts, nwghts, xpins, pins, partvec
       Ouput: partvec, partweights and cut */
    /* Deprecated: use PaToH_Part instead */
    int PaToH_Partition_with_FixCells(PPaToH_Parameters pargs, int _c, int _n,
                                      int *cwghts, int *nwghts, int *xpins,
                                      int *pins, int *partvec, int *partweights,
                                      int *cut);

    /* Input: pargs, _c, _n, _nconst, cwghts, xpins, pins
       Ouput: partvec, partweights and cut */
    /* Deprecated: use PaToH_Part instead */
    int PaToH_MultiConst_Partition(PPaToH_Parameters pargs, int _c, int _n,
                                   int _nconst, int *cwghts,
                                   int *xpins, int *pins, int *partvec,
                                   int *partweights, int *cut);

    /* Unified interface for PaToH with target part weights

       _nconst > 1: calls PaToH_MultiConst_Partition (net weights discarded)
       otherwise    if useFixCells=1  calls PaToH_Partition_with_FixCells
       otherwise calls PaTOH_Partition

       Input: pargs, _c, _n, _nconst, cwghts, nwghts, xpins, pins, targetweights,
       Ouput: partvec, partweights and cut

       if targetweights is NULL each targetweights[i] will be assigned to
       1/pargs->_k

    */
    int PaToH_Part(PPaToH_Parameters pargs, int _c, int _n, int _nconst, int useFixCells,
                   int *cwghts, int *nwghts, int *xpins, int *pins, float *targetweights,
                   int *partvec, int *partweights, int *cut);

/* --- refine --- */
    /* Input: pargs, _c, _n, cwghts, nwghts, xpins, pins, partvec
       Ouput: partvec, partweights and cut */
    int PaToH_Refine_Bisec(PPaToH_Parameters pargs, int _c, int _n,
                           int *cwghts, int *nwghts, int *xpins,
                           int *pins, int *partvec, int *partweights,
                           int *cut);


/* --- utility --- */
    char *PaToH_VersionStr(void);
    int PaToH_Print_Parameter_Abrv(void);

    void PaToH_Print_Parameters(PPaToH_Parameters p);

    int PaToH_Check_User_Parameters(PPaToH_Parameters pargs, int verbose);

    int PaToH_Read_Hypergraph(char *filename, int *_c, int *_n, int *_nconst,
                              int **cwghts, int **nwghts, int **xpins, int **pins);
    int PaToH_Write_Hypergraph(char *filename, int numbering, int _c, int _n, int _nconst,
                               int *cwghts, int *nwghts, int *xpins, int *pins);

    int PaToH_Check_Hypergraph(int _c, int _n, int _nconst,
                               int *cwghts, int *nwghts, int *xpins, int *pins);

    int PaToH_Compute_Cut(int _k, int cuttype, int _c, int _n, int *nwghts,
                          int *xpins, int *pins, int *partvec);

    int PaToH_Compute_Part_Weights(int _k, int _c, int _nconst,
                                   int *cwghts, int *partvec, int *partweights);
    int PaToH_Compute_Part_NetWeights(int _k, int _n, int *nwghts,
                                      int *xpins, int *pins,
                                      int *partvec, int *partinweights);


#ifdef __cplusplus
}
#endif

#endif

