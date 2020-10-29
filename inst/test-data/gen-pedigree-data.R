.do_sim <- TRUE

sim_pedigree_data <- local(envir = new.env(), {
#####
# util function. Do not call these functions. Skip to the end of the file
# and see the plots
invisible(list2env(local({
  idcount <- 0L
  .get_id <- function(n)
    if(n > 0L)
      replicate(n, idcount <<- idcount + 1L) else NULL
  .reset_id <- function()
    (idcount <<- 0L)

  list(.get_id = .get_id, .reset_id = .reset_id)
}), environment()))

.sim_mating <- function(rchild, rmatch, max_depth = 1L, lvl, dadid,
                        momid){
  if(is.null(dadid) || is.null(momid)){
    id <- .get_id(2L)
    father <- mother <- rep(NA_integer_, 2L)
    sex <- 1:2 # 1: male, 2: female
    dadid <- id[1]
    momid <- id[2]
    obslvl <- rep(0L, 2L)
    do_match <- needs_match <- rep(FALSE, 2L)

  } else
    id <- father <- mother <- obslvl <- do_match <- sex <-  NULL

  # sample first lvl
  n_child     <- rchild(1L)
  sex         <- c(sex, (runif(n_child) > .5) + 1L)

  id     <- c(id, .get_id(n_child))
  father <- c(father, rep(dadid, n_child))
  mother <- c(mother, rep(momid, n_child))
  obslvl <- c(obslvl, rep(lvl, n_child))

  do_match    <- c(do_match, rmatch(n_child))
  needs_match <- do_match

  list(id = id, father = father, mother = mother, obslvl = obslvl,
       do_match = do_match, needs_match = needs_match, sex = sex)
}

.sim_household <- function(rchild, rmatch, max_depth, lvl, dadid, momid,
                           max_members = 100L, n_members = 0L){
  out <- list(.sim_mating(rchild, rmatch, max_depth, lvl, dadid, momid))

  get_n_members <- function()
    sum(sapply(out, function(x) length(x[["id"]]))) + n_members
  get_needs_match <- function(obj)
    obj$needs_match & max_depth > obj$obslvl
  shuffle <- function(x)
    if(length(x) == 1L) x else sample(x)

  finish_next <- FALSE
  while(any(
    unlist(sapply(out, get_needs_match))) && !finish_next){
    for(i in shuffle(seq_along(out))){
      if(finish_next)
        break

      for(j in which(get_needs_match(out[[i]]))){
        if(finish_next)
          break

        while(out[[i]]$needs_match[j]){
          # find a matching family
          new_hou <- .sim_household(
            rchild, rmatch, max_depth = out[[i]]$obslvl[j], lvl,
            dadid = NULL, momid = NULL, n_members = get_n_members())
          is_match <- which(
            new_hou$needs_match & new_hou$sex != out[[i]]$sex[j] &
              new_hou$obslvl == out[[i]]$obslvl[j])

          if(length(is_match) == 0L)
            next

          # create the new family
          is_match <- is_match[1L]
          new_hou$needs_match[is_match] <- out[[i]]$needs_match[j] <- FALSE

          if(get_n_members() + length(new_hou$id) >= max_members)
            finish_next <- TRUE # we do not add anymore

          new_pai <- .sim_mating(
            rchild, rmatch, max_depth, lvl = out[[i]]$obslvl[j] + 1L,
            dadid =
              if(out[[i]]$sex[j] == 1L)
                out[[i]]$id[j] else new_hou$id[is_match],
            momid =
              if(out[[i]]$sex[j] == 1L)
                new_hou$id[is_match] else out[[i]]$id[j])

          if(length(new_pai$id) > 0L)
            out <- do.call(c, list(out, list(new_hou), list(new_pai)))
        }

      }
    }
  }

  nam <- names(out[[1L]])
  structure(
    lapply(nam, function(name) do.call(c, lapply(out, "[[", name))),
    names = nam)
}

# simulates a family starting with the oldest members.
#
# Args:
#   rchild: simulation function to sample number of children. Takes an
#           integer with the number of samples to draw.
#   rmatch: simulation function to sample whether a given person meets a
#           partner. Takes an integer with the number of samples to draw.
#   max_depth: max depth of the families.
#   max_members: roughly the amount of members to return at most.
#
# TODO: starting from the bottom might be smarter?
sim_fam <- function(rchild, rmatch, max_depth = 2L, max_members = 100L){
  .reset_id()

  for(i in 1:100){
    out <- as.data.frame(.sim_household(
      rchild, rmatch, max_depth, lvl = 1L, dadid = NULL, momid = NULL,
      max_members = max_members))

    if(NROW(out) > 2L)
      # drop the boring
      break
  }
  old_ids <- out$id
  within(out, {
    id     <- match(id    , old_ids, NA_integer_)
    father <- match(father, old_ids, NA_integer_)
    mother <- match(mother, old_ids, NA_integer_)
    needs_match <- NULL
    do_match <- NULL
  })
}

#####
# parameters. Start with family settings
function(max_depth = 2L, max_members = 100L, sds = 1,
         n_families = 100L, do_plot = FALSE){
  # then the baseline survival function
  base_haz_func <- function(x){
    x <- log(x)
    cbind(x^3, x^2, x)
  }
  d_base_haz_func <- function(x){
    y <- log(x)
    cbind(3 * y^2, 2 * y, 1) / x
  }

  # for a x^3 + b x^2 + c x to monotonically increasing we must have
  # a >= 0 and 2^2 b^2 < 3 * 4 a c <=>  b^2 / 3 a < c
  a_val <- 1e-2
  b_val <- 2e-2
  c_val <- 25 * b_val * b_val / a_val / 3
  omega <- c(a_val, b_val, c_val)
  intecept <- -2.5

  plot(function(x) base_haz_func(exp(x)) %*% omega + intecept,
       xlim = c(log(1e-4), log(10)),
       ylab = expression(-Phi^-1*(S)), xlab = "log(time)")
  plot(function(x) base_haz_func(x) %*% omega + intecept,
       xlim = c(1e-4, 10), ylab = expression(-Phi^-1*(S)), xlab = "time")
  plot(function(x) pnorm(-base_haz_func(x) %*% omega - intecept),
       xlim = c(1e-4, 10), ylab = "S", xlab = "time")
  plot(function(x){
    eta <- drop(-base_haz_func(x) %*% omega - intecept)
    d_eta <- drop(-d_base_haz_func(x) %*% omega)
    -exp(dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE)) * d_eta
  }, xlim = c(1e-4, 10), ylab = "hazard", xlab = "time", ylim = c(0, .02),
  yaxs="i", bty = "l")

  # linear predictor
  gen_x <- function(n)
    cbind(1, rnorm(n))
  beta <- c(intecept, .25)

  # censoring function
  gen_cens <- function(n)
    runif(n, 0, 10)

  # generate families
  library(kinship2)
  dat <- replicate(n_families, {
    for(ii in 1:100){
      dat <- sim_fam(
        rchild = function(n)
          sample.int(size = n, 4L, prob = c(.2, .4, .3, .1)),
        rmatch = function(n) runif(n) > .1,
        max_depth = max_depth, max_members = max_members)
      pedAll <- pedigree(id = dat$id, dadid = dat$father, momid = dat$mother,
                         sex = dat$sex, famid = rep(1, NROW(dat)))["1"]

      rel_mat_full <- t(2  * kinship(pedAll))
      dimnames(rel_mat_full) <- list(dat$id, dat$id)
      rel_mat <- rel_mat_full
      keep <- dat$obslvl == max(dat$obslvl)
      if(sum(keep) > max_members / 10L)
        break
    }

    # matrix for the genetic effect
    rel_mat <- rel_mat[keep, keep, drop = FALSE]
    n_obs <- NROW(rel_mat)

    # simulate the error term
    Sig <- sds[1]^2 * rel_mat
    epsilon <- drop(rnorm(n_obs) %*% chol(Sig))

    Us <- runif(n_obs)
    Zs <- gen_x(n_obs)

    # -Phi^-(S) = b + eta + eps
    # <=> Phi^-(S) + eta + eps = -b
    targets <- qnorm(Us) + drop(Zs %*% beta) + epsilon
    cens <- gen_cens(n_obs)

    # find survival times
    Ys <- mapply(function(x, C){
      f <- function(v)
        drop(-base_haz_func(v) %*% omega) - x
      f_U <- f(C)

      if(f_U > 0)
        return(C * 1.001)

      out <- uniroot(f, interval = c(1e-16, C), f.upper = f_U,
                     tol = .Machine$double.eps^(3/4))
      out$root
    }, x = targets, C = cens)

    is_observed <- Ys < cens
    obs_time <- pmin(Ys, cens)

    list(y = obs_time, event = is_observed, Z = Zs, rel_mat = rel_mat,
         rel_mat_full = rel_mat_full, pedAll = pedAll)
  }, simplify = FALSE)

  # add the true parameter values and save
  list(omega = omega, beta = beta, sds = sds, sim_data = dat)
}
})

if(.do_sim){
  .old_seed <- .Random.seed
  set.seed(2)
  dat <- sim_pedigree_data(do_plot = TRUE, max_members = 20L,
                           n_families = 1000L)
  saveRDS(dat, file.path("inst", "test-data", "pedigree.RDS"))
  .Random.seed <- .old_seed
  rm(.old_seed)

}
