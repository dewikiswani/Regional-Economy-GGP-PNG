library(writexl)

file_path <- "IFPRI_SAM_PNG_2023_SAM.csv" 
output_excel_path <- "PNG_2023_IO_Analysis.xlsx"

# Read CSV
sam_raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)

# Clean account codes and matrix
account_codes <- trimws(as.character(sam_raw$Code))
valid_idx <- !is.na(account_codes) & account_codes != "" & tolower(account_codes) != "total"

account_codes <- account_codes[valid_idx]
sam_mat <- as.matrix(sam_raw[valid_idx, colnames(sam_raw) %in% sam_raw$Code[valid_idx]])

sam_mat <- apply(sam_mat, 2, function(x) as.numeric(gsub(",", "", x)))
sam_mat[is.na(sam_mat)] <- 0
rownames(sam_mat) <- account_codes
colnames(sam_mat) <- account_codes

# Define accounts
act_codes <- account_codes[startsWith(account_codes, "a")]
com_codes <- account_codes[startsWith(account_codes, "c")]
fac_codes <- account_codes[startsWith(account_codes, "f") & account_codes != "ftax"]
tax_codes <- account_codes[account_codes %in% c("ftax", "atax", "mtax", "dtax", "itax")]
hhd_codes <- account_codes[startsWith(account_codes, "hhd")]
gov_code  <- account_codes[account_codes == "gov"]
inv_codes <- account_codes[account_codes %in% c("s-i", "dstk")]
row_code  <- account_codes[account_codes == "row"]

# Make Matrix (V) and Use Matrix (U)
V <- sam_mat[act_codes, com_codes, drop = FALSE]

# Filter active commodities and activities
dom_com_idx   <- which(colSums(V) > 0)
dom_com_codes <- names(dom_com_idx)

active_act_idx <- which(rowSums(V[, dom_com_codes, drop = FALSE]) > 0)
active_act     <- names(active_act_idx)

V_sub <- V[active_act, dom_com_codes, drop = FALSE]
U_sub <- sam_mat[dom_com_codes, active_act, drop = FALSE]

# Control Totals
g <- rowSums(V_sub)               # Gross activity output
q_dom <- colSums(V_sub)           # Domestic commodity output

# 1. Intermediate consumption coefficients: B = U * inv(diag(g))
B <- U_sub %*% diag(1 / ifelse(g == 0, 1, g))
colnames(B) <- active_act
rownames(B) <- dom_com_codes

# 2. Market share transformation: D = V' * inv(diag(g))
# Industry Technology Assumption (ITA): A_c = B %*% D
D <- t(V_sub) %*% diag(1 / ifelse(g == 0, 1, g))
rownames(D) <- dom_com_codes
colnames(D) <- active_act

# Symmetrical Commodity-by-Commodity matrix
A_c <- B %*% D
rownames(A_c) <- dom_com_codes
colnames(A_c) <- dom_com_codes

# Re-scale to intermediate flows Z_c
Z_c <- A_c %*% diag(q_dom)
rownames(Z_c) <- dom_com_codes
colnames(Z_c) <- dom_com_codes

# Factor & Activity Tax allocations
VA_mat <- sam_mat[c(fac_codes, tax_codes), active_act, drop = FALSE]
VA_c   <- VA_mat %*% diag(1 / ifelse(g == 0, 1, g)) %*% D %*% diag(q_dom)
rownames(VA_c) <- rownames(VA_mat)
colnames(VA_c) <- dom_com_codes

# Final Demand Categories
FD_hhd <- rowSums(sam_mat[dom_com_codes, hhd_codes, drop = FALSE])
FD_gov <- rowSums(sam_mat[dom_com_codes, gov_code, drop = FALSE])
FD_inv <- rowSums(sam_mat[dom_com_codes, inv_codes, drop = FALSE])
FD_exp <- sam_mat[dom_com_codes, row_code]
M_c    <- sam_mat[row_code, dom_com_codes]

Final_Demand <- cbind(
  Households = FD_hhd,
  Government = FD_gov,
  Investment = FD_inv,
  Exports    = FD_exp
)

Total_Gross_Output <- rowSums(Z_c) + rowSums(Final_Demand) - as.vector(M_c)
Total_Gross_Input  <- colSums(Z_c) + colSums(VA_c)

# Construct Table
top_block <- cbind(Z_c, Final_Demand, Imports = -as.vector(M_c), Total_Output = Total_Gross_Output)
bottom_block <- cbind(VA_c, matrix(0, nrow = nrow(VA_c), ncol = ncol(Final_Demand) + 2))
colnames(bottom_block) <- colnames(top_block)

io_table <- rbind(top_block, bottom_block)
io_table <- rbind(io_table, Total_Input = c(Total_Gross_Input, rep(NA, ncol(Final_Demand) + 2)))

# Leontief Inverse: L = (I - A_c)^-1
I_mat <- diag(nrow(A_c))
L_c   <- solve(I_mat - A_c)
rownames(L_c) <- dom_com_codes
colnames(L_c) <- dom_com_codes

#Finalized the data and construct I-O Table,incorporating Import

io_table_df <- data.frame(Account = rownames(io_table), round(io_table, 6), check.names = FALSE)

col_at_idx <- which(names(io_table_df) %in% c("AT", "Imports"))[1]
if (is.na(col_at_idx)) col_at_idx <- 46

col_au_idx <- which(names(io_table_df) %in% c("AU", "Total_Output"))[1]
if (is.na(col_au_idx)) col_au_idx <- 47

row_50_idx <- which(io_table_df[[1]] %in% c("Total_Input", "Total Input"))[1]
if (is.na(row_50_idx)) row_50_idx <- 49

# Number of production sectors/commodities (rows 1 to 40)
n_sectors <- which(io_table_df[[1]] == "cosrv")[1]
if (is.na(n_sectors)) n_sectors <- 40

# -------------------------------------------------------------------------
# 2. Compute Current Sums
# -------------------------------------------------------------------------
sum_col_au <- sum(as.numeric(io_table_df[1:n_sectors, col_au_idx]), na.rm = TRUE)
sum_row_50 <- sum(as.numeric(io_table_df[row_50_idx, 2:(n_sectors + 1)]), na.rm = TRUE)

cat("Sum of Column AU (Total Output):", sum_col_au, "\n")
cat("Sum of Row 50 (Total Input)    :", sum_row_50, "\n")

diff_val <- sum_row_50 - sum_col_au
cat("Discrepancy (Row 50 - Col AU)  :", diff_val, "\n\n")

if (abs(diff_val) > 1e-6) {
  at_values <- as.numeric(io_table_df[1:n_sectors, col_at_idx])
  sum_at <- sum(at_values, na.rm = TRUE)
  
  # Allocate difference proportionally based on the magnitude of each row's imports
  # New AT_sum = sum_at + diff_val
  scaling_ratio <- (sum_at + diff_val) / sum_at
  
  new_at_values <- at_values * scaling_ratio
  
  # Ensure values remain non-positive
  if (any(new_at_values > 0, na.rm = TRUE)) {
    warning("Some values in Column AT became positive. Please inspect allocation weights.")
  }
  
  # Update column AT
  io_table_df[1:n_sectors, col_at_idx] <- new_at_values
  
  # Recalculate Column AU (Total_Output = Sum of Intermediate & Final Demand + Imports)
  # Columns 2 to AT (index col_at_idx)
  io_table_df[1:n_sectors, col_au_idx] <- rowSums(
    sapply(io_table_df[1:n_sectors, 2:col_at_idx], as.numeric), 
    na.rm = TRUE
  )
  
  # Verify new sums
  updated_sum_au <- sum(as.numeric(io_table_df[1:n_sectors, col_au_idx]), na.rm = TRUE)
  cat("=== After Reallocation ===\n")
  cat("Updated Sum of Column AT:", sum(new_at_values, na.rm = TRUE), "\n")
  cat("Updated Sum of Column AU:", updated_sum_au, "\n")
  cat("Sum of Row 50           :", sum_row_50, "\n")
  cat("Difference              :", sum_row_50 - updated_sum_au, "\n")
} else {
  cat("Column AU and Row 50 are already balanced.\n")
}

# Export to Excel for PNG
io_table_export <- io_table_df

# 2. Convert all columns except 'Account' to numeric and round to 6 decimal digits
numeric_cols <- names(io_table_export)[-1]  # all columns except the first ('Account')

io_table_export[numeric_cols] <- lapply(io_table_export[numeric_cols], function(col) {
  num_col <- suppressWarnings(as.numeric(col))
  # Replace NA with 0 if any blanks exist
  num_col[is.na(num_col)] <- 0
  round(num_col, 6)
})

write_xlsx(
  x = list("PNG_2023_IO" = io_table_export),
  path = "PNG_2023_IO_Table.xlsx",
  col_names = TRUE,
  format_headers = TRUE
)

cat("Successfully exported 'PNG_2023_IO_Table.xlsx' with 6 decimal places.\n")

# ==============================================================================
# Regionalization of PNG 2019 I-O Table to Oro 
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup Sector Structure and Baseline Data
# ------------------------------------------------------------------------------
n_sec <- 40  # 40 production sectors (cmaiz to cosrv)
sec_names <- as.character(io_table_df[1:n_sec, 1])

# National Intermediate Transaction Matrix Z^N (40x40)
Z_nat <- as.matrix(sapply(io_table_df[1:n_sec, 2:(n_sec + 1)], as.numeric))
colnames(Z_nat) <- sec_names
rownames(Z_nat) <- sec_names

# National Value Added & Taxes Matrix V^N (8 rows x 40 cols)
# Rows 41:48 (flab-n, flab-p, flab-s, flnd, fcap, dtax, ftax, mtax)
V_nat <- as.matrix(sapply(io_table_df[(n_sec + 1):(n_sec + 8), 2:(n_sec + 1)], as.numeric))
rownames(V_nat) <- as.character(io_table_df[(n_sec + 1):(n_sec + 8), 1])
colnames(V_nat) <- sec_names

# National Final Demand (Households, Govt, Investment, Exports) (40x4)
col_fd_start <- which(names(io_table_df) == "Households")[1]
if (is.na(col_fd_start)) col_fd_start <- 42  # Column AP in 1-based index

FD_nat <- as.matrix(sapply(io_table_df[1:n_sec, col_fd_start:(col_fd_start + 3)], as.numeric))
colnames(FD_nat) <- c("Households", "Government", "Investment", "Exports")

# National Total Output X^N (Column AU / index 47)
col_au_idx <- which(names(io_table_df) %in% c("AU", "Total_Output"))[1]
if (is.na(col_au_idx)) col_au_idx <- 47
X_nat <- as.numeric(io_table_df[1:n_sec, col_au_idx])
names(X_nat) <- sec_names

X_nat_safe <- ifelse(X_nat == 0, 1e-6, X_nat)

# ------------------------------------------------------------------------------
# 2. Empirical Sectoral Output Shares for Oro Province
# ------------------------------------------------------------------------------
# Baseline population share from NSO 2011 / 2024 Census: ~2.4%
pop_share <- 0.024
reg_output_share <- setNames(rep(pop_share, n_sec), sec_names)

# Calibrate specific sectors based on commodity and administrative records:
# WARNING: ALL OF THIS ASSUMPTION WILL REQUIRES VALIDATION!
# A. Agriculture & Primary Commodities
reg_output_share["coils"] <- 0.150  # Oil palm (OPIC / NBPOL Higaturu Estate)
reg_output_share["cfore"] <- 0.045  # Forestry (PNGFA active concessions)
reg_output_share["cfish"] <- 0.035  # Marine / coastal fisheries
reg_output_share["ccoff"] <- 0.030  # Coffee (CIC - Managalas/Kokoda smallholders)
reg_output_share["cocrp"] <- 0.040  # Cocoa & rubber (Cocoa Board estimates)
reg_output_share["csugr"] <- 0.005  # Ramu Sugar is in Madang/Morobe

# B. Mining & Extraction
reg_output_share["cmine"] <- 0.001  # No commercial active mine in Oro

# C. Heavy Manufacturing & Refining (almost entirely imported from Lae/POM)
reg_output_share["cmach"] <- 0.002
reg_output_share["cmetl"] <- 0.002
reg_output_share["cchem"] <- 0.003
reg_output_share["ctext"] <- 0.005
reg_output_share["cnmet"] <- 0.005

# Compute Regional Output Vector (X^R)
X_reg <- X_nat * reg_output_share
tot_X_nat <- sum(X_nat)
tot_X_reg <- sum(X_reg)
R_size <- tot_X_reg / tot_X_nat

# ------------------------------------------------------------------------------
# 3. Flegg's Location Quotient (FLQ) Regionalization
# ------------------------------------------------------------------------------
# Simple Location Quotient (SLQ)
SLQ <- reg_output_share / R_size

# Flegg parameter (delta = 0.3 for provincial scale)
delta <- 0.3
lambda_star <- (log2(1 + R_size))^delta
FLQ <- SLQ * lambda_star

# Retention factor q_i (bounded at 1.0)
q_factor <- pmin(FLQ, 1.0)

# National Technical Coefficients A^N
A_nat <- sweep(Z_nat, 2, X_nat_safe, FUN = "/")

# Intra-regional Technical Coefficients A^R
A_reg_initial <- sweep(A_nat, 1, q_factor, FUN = "*")

# Initial Regional Intermediate Matrix Z^R
Z_reg_initial <- sweep(A_reg_initial, 2, X_reg, FUN = "*")

# ------------------------------------------------------------------------------
# 4. Regional Value Added and Final Demand Setup
# ------------------------------------------------------------------------------
# Regional Value Added (proportional to sector output)
V_reg <- sweep(sweep(V_nat, 2, X_nat_safe, FUN = "/"), 2, X_reg, FUN = "*")

# Regional Final Demand (Households scale with pop_share; Exports scale with sectoral output)
FD_reg <- FD_nat
FD_reg[, c("Households", "Government", "Investment")] <- FD_nat[, c("Households", "Government", "Investment")] * pop_share
FD_reg[, "Exports"] <- FD_nat[, "Exports"] * reg_output_share

# ------------------------------------------------------------------------------
# 5. Target Margins & RAS Biproportional Balancing
# ------------------------------------------------------------------------------
# Column target for intermediate inputs:
# u_target_j = X_reg_j - sum(V_reg[, j])
u_target <- pmax(X_reg - colSums(V_reg), 0)

# Row target for intermediate output:
# r_target_i = X_reg_i - sum(FD_reg[i, ])
# Any gap where local final demand exceeds output is handled via external imports
net_fd <- rowSums(FD_reg)
r_target <- pmax(X_reg - net_fd, 0)

# RAS Function
ras_balance <- function(Z, target_r, target_u, max_iter = 1000, tol = 1e-6) {
  Z_k <- Z
  for (iter in 1:max_iter) {
    # 1. Row scaling (r step)
    current_row <- rowSums(Z_k)
    current_row_safe <- ifelse(current_row == 0, 1e-6, current_row)
    r_mult <- ifelse(target_r == 0, 0, target_r / current_row_safe)
    Z_k <- sweep(Z_k, 1, r_mult, FUN = "*")
    
    # 2. Column scaling (s step)
    current_col <- colSums(Z_k)
    current_col_safe <- ifelse(current_col == 0, 1e-6, current_col)
    s_mult <- ifelse(target_u == 0, 0, target_u / current_col_safe)
    Z_k <- sweep(Z_k, 2, s_mult, FUN = "*")
    
    # Check convergence
    row_err <- max(abs(rowSums(Z_k) - target_r))
    col_err <- max(abs(colSums(Z_k) - target_u))
    
    if (max(row_err, col_err) < tol) {
      cat(sprintf("RAS successfully converged at iteration %d.\n", iter))
      break
    }
  }
  return(Z_k)
}

# Execute RAS
Z_reg_balanced <- ras_balance(Z_reg_initial, r_target, u_target)

# ------------------------------------------------------------------------------
# 6. Corrected Assembly: Preserving Strictly Negative Imports
# ------------------------------------------------------------------------------
# Total demand prior to trade adjustment
tot_demand <- rowSums(Z_reg_balanced) + rowSums(FD_reg)

# Calculate trade balance for each sector
# If balance < 0 (Demand > Output) -> Oro needs IMPORTS (negative value)
# If balance > 0 (Output > Demand) -> Oro has REGIONAL EXPORTS to rest of PNG
trade_balance <- X_reg - tot_demand

# Initialize imports as non-positive
M_reg_corrected <- rep(0, n_sec)
names(M_reg_corrected) <- sec_names

# Deficit sectors: imported from external / rest of PNG (strictly <= 0)
M_reg_corrected[trade_balance < 0] <- trade_balance[trade_balance < 0]

# Surplus sectors: exported to other provinces -> add to the Exports column
FD_reg_corrected <- FD_reg
FD_reg_corrected[trade_balance > 0, "Exports"] <- 
  FD_reg_corrected[trade_balance > 0, "Exports"] + trade_balance[trade_balance > 0]

# Verify row sum identity: Total Output = Intermediate + Final Demand + Imports
recalc_output <- rowSums(Z_reg_balanced) + rowSums(FD_reg_corrected) + M_reg_corrected
cat("Max output calculation discrepancy:", max(abs(recalc_output - X_reg)), "\n")
cat("Any positive values in Imports?    :", any(M_reg_corrected > 0), "\n")

# Assemble the final table
oro_io_table <- data.frame(
  Account = sec_names,
  round(Z_reg_balanced, 4),
  round(FD_reg_corrected, 4),
  Imports = round(M_reg_corrected, 4),
  Total_Output = round(X_reg, 4),
  check.names = FALSE
)

# Append Value Added rows
va_rows <- data.frame(
  Account = rownames(V_reg),
  round(V_reg, 4),
  matrix(0, nrow = nrow(V_reg), ncol = ncol(FD_reg_corrected) + 2, 
         dimnames = list(NULL, c(colnames(FD_reg_corrected), "Imports", "Total_Output"))),
  check.names = FALSE
)
colnames(va_rows) <- colnames(oro_io_table)

# Append Total Input row
total_input_vec <- c(
  "Total_Input",
  round(colSums(Z_reg_balanced) + colSums(V_reg), 4),
  rep(0, ncol(FD_reg_corrected) + 2)
)

oro_io_table <- rbind(oro_io_table, va_rows, total_input_vec)

# ------------------------------------------------------------------------------
# 7. Verification of Balance
# ------------------------------------------------------------------------------
chk_output <- as.numeric(oro_io_table[1:n_sec, "Total_Output"])
chk_input  <- as.numeric(oro_io_table[nrow(oro_io_table), 2:(n_sec + 1)])
max_imbalance <- max(abs(chk_output - chk_input))

cat("=== Final Validation ===\n")
cat("Total Regional Output (Oro):", sum(chk_output), "\n")
cat("Total Regional Input  (Oro):", sum(chk_input), "\n")
cat("Maximum Sectoral Discrepancy:", max_imbalance, "\n")


# -------------------------------------------------------------------------
# FINALIZATIOn I-O
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# 1. Identify Target Dimensions and Sector Rows
# -------------------------------------------------------------------------
col_at_idx <- which(names(oro_io_table) %in% c("AT", "Imports"))[1]
if (is.na(col_at_idx)) col_at_idx <- 46

col_au_idx <- which(names(oro_io_table) %in% c("AU", "Total_Output"))[1]
if (is.na(col_au_idx)) col_au_idx <- 47

row_50_idx <- which(oro_io_table[[1]] %in% c("Total_Input", "Total Input"))[1]
if (is.na(row_50_idx)) row_50_idx <- 49

# Number of production sectors/commodities (rows 1 to 40)
n_sectors <- which(oro_io_table[[1]] == "cosrv")[1]
if (is.na(n_sectors)) n_sectors <- 40

# -------------------------------------------------------------------------
# 2. Compute Current Sums
# -------------------------------------------------------------------------
sum_col_au <- sum(as.numeric(oro_io_table[1:n_sectors, col_au_idx]), na.rm = TRUE)
sum_row_50 <- sum(as.numeric(oro_io_table[row_50_idx, 2:(n_sectors + 1)]), na.rm = TRUE)

cat("Sum of Column AU (Total Output):", sum_col_au, "\n")
cat("Sum of Row 50 (Total Input)    :", sum_row_50, "\n")

diff_val <- sum_row_50 - sum_col_au
cat("Discrepancy (Row 50 - Col AU)  :", diff_val, "\n\n")

if (abs(diff_val) > 1e-6) {
  at_values <- as.numeric(oro_io_table[1:n_sectors, col_at_idx])
  sum_at <- sum(at_values, na.rm = TRUE)
  
  # Allocate difference proportionally based on the magnitude of each row's imports
  # New AT_sum = sum_at + diff_val
  scaling_ratio <- (sum_at + diff_val) / sum_at
  
  new_at_values <- at_values * scaling_ratio
  
  # Ensure values remain non-positive
  if (any(new_at_values > 0, na.rm = TRUE)) {
    warning("Some values in Column AT became positive. Please inspect allocation weights.")
  }
  
  # Update column AT
  oro_io_table[1:n_sectors, col_at_idx] <- new_at_values
  
  # Recalculate Column AU (Total_Output = Sum of Intermediate & Final Demand + Imports)
  # Columns 2 to AT (index col_at_idx)
  oro_io_table[1:n_sectors, col_au_idx] <- rowSums(
    sapply(oro_io_table[1:n_sectors, 2:col_at_idx], as.numeric), 
    na.rm = TRUE
  )
  
  # Verify new sums
  updated_sum_au <- sum(as.numeric(oro_io_table[1:n_sectors, col_au_idx]), na.rm = TRUE)
  cat("=== After Reallocation ===\n")
  cat("Updated Sum of Column AT:", sum(new_at_values, na.rm = TRUE), "\n")
  cat("Updated Sum of Column AU:", updated_sum_au, "\n")
  cat("Sum of Row 50           :", sum_row_50, "\n")
  cat("Difference              :", sum_row_50 - updated_sum_au, "\n")
} else {
  cat("Column AU and Row 50 are already balanced.\n")
}

# Export to Excel
oro_io_export <- oro_io_table

# 2. Convert all columns except 'Account' to numeric and round to 6 decimal digits
numeric_cols <- names(oro_io_export)[-1]  # all columns except the first ('Account')

oro_io_export[numeric_cols] <- lapply(oro_io_export[numeric_cols], function(col) {
  num_col <- suppressWarnings(as.numeric(col))
  # Replace NA with 0 if any blanks exist
  num_col[is.na(num_col)] <- 0
  round(num_col, 6)
})

write_xlsx(
  x = list("Oro_2023_IO" = oro_io_export),
  path = "Oro_2023_IO_Table.xlsx",
  col_names = TRUE,
  format_headers = TRUE
)

cat("Successfully exported 'Oro_2023_IO_Table.xlsx' with 6 decimal places.\n")



