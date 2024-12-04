#### CODE TO GENERATE JOBS FOR CLUSTER ####

# Settings ----------------------------------------------------------------
subcode <- NA

for (g in seq(1, 100)) {

  mycode <- paste(
    "#!/bin/bash\n",
    "#SBATCH -n 1\n",
    "#SBATCH -N 1\n",
    "#SBATCH -t 0-00:02:00\n",
    "#SBATCH -p shared\n",
    "#SBATCH --mem=10G\n",
    "#SBATCH -c 15\n",
    "#SBATCH -J Jobs/job_", g, "\n",
    "#SBATCH -o Jobs/job_", g, ".out\n",
    "#SBATCH -e Jobs/job_", g, ".err\n",
    # "#SBATCH --test-only\n", # sanity check at beginning
    "\n",
    "module load R\n",
    "module load gcc/12.2.0-fasrc01\n",
    "unset R_LIBS_SITE\n",
    "export R_LIBS_USER=$HOME/apps/R/4.2.2:$R_LIBS_USER\n",
    "Rscript run_log_reg_simulation.R", " ",
    g, " ",
    sep = ""
  )

  myname <- paste("Jobs/job_",
                  g,
                  ".sh",
                  sep="")

  subcode <- paste(subcode,
                   "sbatch ",
                   myname,
                   ".sh\n",
                   sep="")

  write.table(mycode,
              myname,
              quote=F,
              row.names = F,
              col.names = F)

}



write.table(subcode,
            "Jobs/_subcode.txt",
            quote = F,
            row.names = F,
            col.names = F)

