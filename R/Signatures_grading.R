#' @title Giving immunotherapy prognosis using MEMTS.
#' @description The function will return a vector calculated by Metastasis Related Epithelial-Mesenchymal Transition Signature.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

MEMTS_grading <- function(exp_mtr){
  MEMTS <- c('SNAI2','PFN2','NOTCH2','NID2','MEST','MATN2','LAMA1','ITGB3','GPX7','FBN2','ECM2','DPYSL3','BDNF')
  Expr <- dataPreprocess(exp_mtr,MEMTS,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Giving immunotherapy prognosis using PRGScore.
#' @description The function will return a vector calculated using pyroptosis-related gene score.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

PRGScore_grading <- function(exp_mtr){
  PRGScore <- c('GSDMB','TNFRSF14','OAS1','HSH2D','APOL2','APOL1','HCG4P11','RP1-50J22.4','CD96','RP1-102E24.10','APOBEC3G','CYP4F12','ETV7','FBXO6','B3GNT3','APOBEC3D','HLA-K','OR2I1P','EXOC3L4','PSMB10','RP11-291B21.2','HCG4P7','ZNF683','GPR114','TJP3','UBD','LGALS9','CGREF1','TMC4','FOXA1','TRGC2','RP11-876N24.3','INPP1','AC092580.4','NRIR','UBA7','PSMB8','APOL6','CCDC64B','KIR2DL4','FXYD3','HLA-G','PEG10','USP30-AS1','PSMB8-AS1','TRBV9','AGR2','AC002331.1','TRGV10','ELF3','RTP4','PRR9','HLA-F','PSME2','REC8','CTSS','HLA-E')
  Expr <- dataPreprocess(exp_mtr,PRGScore,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Giving immunotherapy prognosis using Angiogenesis.
#' @description The function will return a vector calculated using Angiogenesis.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

Angiogenesis_grading <- function(exp_mtr){
  Angiogenesis <- c('VEGFA','KDR','ESM1','PECAM1','ANGPTL4','CD34')
  Expr <- dataPreprocess(exp_mtr,Angiogenesis,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Giving immunotherapy prognosis using T-effector.
#' @description The function will return a vector calculated using T-effector.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

Teffector_grading <- function(exp_mtr){
  Teffector <- c('CD8A','EOMES','PRF1','IFNG','CD274')
  Expr <- dataPreprocess(exp_mtr,Teffector,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Giving immunotherapy prognosis using myeloid inflammatory.
#' @description The function will return a vector calculated using myeloid inflammatory.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

Myeloid_inflammatory_grading <- function(exp_mtr){
  MI <- c('IL-6','CXCL1','CXCL2','CXCL3','CXCL8','PTGS2')
  Expr <- dataPreprocess(exp_mtr,MI,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}


#' @title Giving immunotherapy prognosis using IFNG.
#' @description The function will return a vector calculated using IFNG.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

IFNG_grading <- function(exp_mtr){
  IFNG <- c('CTLA4','STAT1','IRF1','TAP2','GBP2','HLA-DRB5')
  Expr <- dataPreprocess(exp_mtr,IFNG,turn2HL = FALSE)
  result <- apply(Expr, 2, mean)
  names(result) <- colnames(Expr)
  return(result)
}

#' @title Giving immunotherapy prognosis using DDR.
#' @description The function will return a vector calculated using DDR.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

DDR_grading <- function(exp_mtr){
  DDR <- c('UNG','SMUG1','MBD4','OGG1','MUTYH','NTHL1','MPG','NEIL1','NEIL2','NEIL3','APEX1','APEX2','LIG3','XRCC1','PNKP','APLF','PARP1','PARP2','PARP3','MGMT','ALKBH2','ALKBH3','TDP1','TDP2','MSH2','MSH3','MSH6','MLH1','PMS2','MSH4','MSH5','MLH3','PMS1','XPC','RAD23B','CETN2','RAD23A','XPA','DDB1','DDB2','RPA1','RPA2','RPA3','ERCC3','ERCC2','GTF2H1','GTF2H2','GTF2H3','GTF2H4','GTF2H5','CDK7','CCNH','MT1','ERCC5','ERCC1','ERCC4','LIG1','ERCC8','ERCC6','UVSSA','XAB2','MMS19','RAD51','RAD51B','RAD51D','DMC1','XRCC2','XRCC3','RAD52','RAD54L','RAD54B','BRCA1','SHFM1','RAD50','MRE11A','NBN','RBBP8','MUS81','EME1','EME2','GEN1','FANCA','FANCB','FANCC','BRCA2','FANCD2','FANCE','FANCF','FANCG','FANCI','BRIP1','FANCL','FANCM','PALB2','RAD51C','XRCC6','XRCC5','PRKDC','LIG4','XRCC4','DCLRE1C','NHEJ1','NUDT1','DUT','RRM2B','POLB','POLG','POLD1','POLE','PC','REV3L','MAD2L2','POLH','POLI','POLQ','POLK','POLL','POLM','POLN','FEN1','FAN1','TREX1','EXO1','APTX','ENDOV','UBE2A','UBE2B','RAD18','SHPRH','HLTF','RNF168','SPRTN','RNF8','RNF4','UBE2V2','UBE2N','H2AFX','CHAF1A','SETMAR','BLM','WRN','RECQL4','ATM','DCLRE1A','DCLRE1B','RPA4','PRPF19','RECQL','RECQL5','HELQ','RDM1','ATR','ATRIP','MDC1','RAD1','RAD9A','HUS1','RAD17','CHEK1','CHEK2','TP53','TP53BP1','RIF1','TOPBP1','CLK2','PER1')
  Expr <- dataPreprocess(exp_mtr, DDR, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using CD_8_T_effector.
#' @description The function will return a vector calculated using CD_8_T_effector.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

CD8Teffector_grading <- function(exp_mtr){
  CD_8_T_effector <- c('CD8A','GZMA','GZMB','IFNG','CXCL9','CXCL10','PRF1','TBX21')
  Expr <- dataPreprocess(exp_mtr, CD_8_T_effector, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using CellCycle_Reg.
#' @description The function will return a vector calculated using CellCycle_Reg.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

CellCycleReg_grading <- function(exp_mtr){
  CellCycle_Reg <- c('ATM','CDKN1A','CDKN2A','MDM2','TP53','CCND1','RB1','CCNE1','FBXW7','E2F3')
  Expr <- dataPreprocess(exp_mtr, CellCycle_Reg, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using Pan_F_TBRs.
#' @description The function will return a vector calculated using Pan_F_TBRs.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

PanFTBRs_grading <- function(exp_mtr){
  Pan_F_TBRs <- c('ACTA2','ACTG2','ADAM12','ADAM19','CNN1','COL4A1','CTGF','CTPS1','FAM101B','FSTL3','HSPB1','IGFBP3','PXDC1','SEMA7A','SH3PXD2A','TAGLN','TGFBI','TNS1','TPM1')
  Expr <- dataPreprocess(exp_mtr, Pan_F_TBRs, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using EMT1.
#' @description The function will return a vector calculated using EMT1.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

EMT1_grading <- function(exp_mtr){
  EMT1 <- c('CLDN3','CLDN7','CLDN4','CDH1','VIM','TWIST1','ZEB1','ZEB2')
  Expr <- dataPreprocess(exp_mtr, EMT1, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using EMT2.
#' @description The function will return a vector calculated using EMT2.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

EMT2_grading <- function(exp_mtr){
  EMT2 <- c('AXL','ROR2','WNT5A','LOXL2','TWIST2','TAGLN','FAP')
  Expr <- dataPreprocess(exp_mtr, EMT2, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using EMT3.
#' @description The function will return a vector calculated using EMT3.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

EMT3_grading <- function(exp_mtr){
  EMT3 <- c('SOX9','TWIST1','FOXF1','ZEB1','ZEB2','GATA6')
  Expr <- dataPreprocess(exp_mtr, EMT3, turn2HL = FALSE)
  average <- apply(Expr, 1, mean)
  standard_error <- apply(Expr, 1, sd)
  ZScore <- (Expr - average) / standard_error
  result  <- stats::prcomp(ZScore, center = F, scale = F)$rotation[,1]
  return(result)
}


#' @title Giving immunotherapy prognosis using MSKCC.
#' @description The function will return a vector calculated using MSKCC.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

MSKCC_grading <- function(exp_mtr){
  MSKCC <- c('TP53','PIK3CA','ATM')
  Expr <- dataPreprocess(exp_mtr, MSKCC, turn2HL = FALSE)
  if(nrow(Expr) == 3){
    result  <- -0.492 * Expr[1,] + 0.562 * Expr[2,] - 1.454 * Expr[3,]
  } else{
    warning("There absence in gene TP53, PIK3CA or ATM. Please check your data!")
    return()
    }
  return(result)
}


#' @title Giving immunotherapy prognosis using LMRGPI.
#' @description The function will return a vector calculated using LMRGPI.
#' @param exp_mtr An expression matrix which rownames are gene symbol and colnames are sample ID.
#' @export
#'

LMRGPI_grading <- function(exp_mtr){
  LMRGPI <- c('ANGPTL4','NPAS2','SLCO1B3','ACOXL','ALOX15','B3GALNT1')
  Expr <- dataPreprocess(exp_mtr, LMRGPI, turn2HL = FALSE)
  if(nrow(Expr) == 6){
    result  <- 108.2 * Expr[1,] + 265.1 * Expr[2,] - 83 * Expr[3,] - 261.15 * Expr[4,] -191.3 * Expr[5,] + 177 *  Expr[6,]
  } else{
    warning("There absence in LMRGPI genes. Please check your data!")
    return()
  }
  return(result)
}




