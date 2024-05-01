#' @title dd
#' @param mtr description
#' @export

cell_name_filter <- function(mtr){
  rn <- rownames(mtr)
  rn <- sub("T\\.cell\\.CD8.","T cells CD8",rn) %>%
    sub("T\\.cells\\.CD8","T cells CD8",.) %>%
    sub("T_cells_CD8","T cells CD8",.,fixed = TRUE) %>%
    sub("T cells regulatory \\(Tregs\\)","Tregs",.) %>%
    sub("T\\_regulatory\\_cells","Tregs",.) %>%
    sub("T\\.cell\\.CD4.","T cells CD4",.) %>%
    sub("T\\.cells\\.CD4","T cells CD4",.) %>%
    sub("T\\_cells\\_CD4","T cells CD4",.) %>%
    sub("T cells gamma delta","Yd T cells",.) %>%
    sub("T\\_cells\\_gamma\\_delta","Yd T cells",.) %>%
    sub("^B.cell$","B cells",.) %>%
    sub("B\\-cells","B cells",.) %>%
    sub("B.cells","B cells",.,fixed = TRUE) %>%
    sub("B\\_cells","B cells",.) %>%
    sub("Plasma\\_cells","Plasma cells",.) %>%
    sub("Myeloid\\.dentritic\\.cell","mDCs",.) %>%
    sub("Myeloid dendritic cells","mDCs",.) %>%
    sub("Mast\\_cells","Mast cells",.) %>%
    sub("NK\\_cells","NK cells",.) %>%
    sub("Dendritic cells","DCs",.) %>%
    sub("Dendritic\\.cells","DCs",.) %>%
    sub("Cytotoxic lymphocytes","CTLs",.) %>%
    sub("T cells CD4 memory activated","CD4+ T memory activated",.) %>%
    sub("Immune\\_Score","Immune Score",.)
  rownames(mtr) <- rn
  mtr
}
