# This framework provides flexible parameters to tune for influence maximization problem

library(igraph)
library(dplyr)
library(lubridate)
library(influence.mining)
library(jsonlite)
library(clipr)
library(parallel)
#library(doMC) # For linux
library(doSNOW) # For Windows
library(future.apply)
root_dir <- "imine-framework-evaluation"

load_graphs <- function() {
  karate <- largest_component(read.graph("dataset/karate_club.txt", directed=FALSE))
  karate$name <- "karate"
  football <- largest_component(read.graph("dataset/american_football.txt", directed=FALSE))
  football$name <- "football"
  world_trade <- largest_component(read.graph("dataset/world_trade.txt", directed=FALSE))
  world_trade$name <- "world_trade"
  nematode <- largest_component(read.graph("dataset/nematode_neural_network.txt", directed=FALSE))
  nematode$name <- "nematode"
  automata <- largest_component(read.graph("dataset/automata_nlp.txt", directed=FALSE))
  automata$name <- "automata"
  author <- largest_component(read.graph("dataset/author_netscience.txt", directed=FALSE))
  author$name <- "author"
  protein <- largest_component(read.graph("dataset/protein_barabasi.txt", directed=FALSE))
  protein$name <- "protein"
  influence_citation <- largest_component(read.graph("dataset/influence_citation_network.txt", directed=FALSE, format = "ncol"))
  influence_citation$name <- "influence_citation"
  routers <- largest_component(read.graph("dataset/tech-routers-rf.txt", directed=FALSE, format = "ncol"))
  routers$name <- "tech_routers"
  hamster <- largest_component(read.graph("dataset/petster-hamster-friend.txt", directed=FALSE, format = "ncol"))
  hamster$name <- "hamster"
  moreno <- largest_component(read.graph("dataset/moreno_health.txt", directed=FALSE, format = "ncol"))
  moreno$name <- "moreno"
  fb_tvshows <- largest_component(read.graph("dataset/fb-pages-tvshow-edges.txt", directed=FALSE, format = "ncol"))
  fb_tvshows$name <- "fb_tvshows"
  political <- largest_component(read.graph("dataset/political_blog.txt", directed=FALSE, format = "ncol"))
  political$name <- "political_blog"
  usairport <- largest_component(read.graph("dataset/opsahl_usairport.txt", directed=FALSE, format = "ncol"))
  usairport$name <- "usairport"
  spam <- largest_component(read.graph("dataset/web-spam.txt", directed=FALSE, format = "ncol"))
  spam$name <- "web_spam"
  caida <- largest_component(read.graph("dataset/as-caida.txt", directed=FALSE))
  caida$name <- "caida"
  arxiv <- largest_component(read.graph("dataset/arxiv_collaboration.txt", directed=FALSE, format = "ncol"))
  arxiv$name <- "arxiv"
  ita2000 <- largest_component(read.graph("dataset/ita2000.txt", directed=FALSE))
  ita2000$name <- "ita2000"
  jdk6 <- largest_component(read.graph("dataset/jdk6_dependencies.txt", directed=FALSE))
  jdk6$name <- "jdk6_dependencies"
  johnshopkins <- largest_component(read.graph("dataset/socfb-JohnsHopkins55.txt", directed=FALSE))
  johnshopkins$name <- "fb_johnshopkins"
  linux <- largest_component(read.graph("dataset/linux_code.txt", directed=FALSE))
  linux$name <- "linux_code"
  city2010 <- largest_component(read.graph("dataset/citycity_weighted_2010.csv", directed=FALSE, format = "ncol"))
  city2010$name <- "city2010"
  V(city2010)$name <- 1:vcount(city2010)
  city2013 <- largest_component(read.graph("dataset/citycity_weighted_2013.csv", directed=FALSE, format = "ncol"))
  city2013$name <- "city2013"
  V(city2013)$name <- 1:vcount(city2013)
  city2016 <- largest_component(read.graph("dataset/citycity_weighted_2016.csv", directed=FALSE, format = "ncol"))
  city2016$name <- "city2016"
  V(city2016)$name <- 1:vcount(city2016)
  city2019 <- largest_component(read.graph("dataset/citycity_weighted_2019.csv", directed=FALSE, format = "ncol"))
  city2019$name <- "city2019"
  V(city2019)$name <- 1:vcount(city2019)
  list(karate, football, world_trade, nematode, automata, author, protein, influence_citation, routers, hamster, moreno, fb_tvshows, political, usairport, spam, caida, arxiv, ita2000, jdk6, johnshopkins, linux, youtube, epinions, wordnet, city2010, city2013, city2016, city2019)
}


#' Influence mining framework implementation based on Adaptive centrality method
#'
#' @name adaptive_centrality_influence_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @param prob probability threshold for INFLUENCE_LT. Default is 0.1
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
adaptive_centrality_influence_framework <- function(graph, criteria=c("TIME","BUDGET"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), centrality_method=c("DEGREE", "BETWEENNESS", "CLOSENESS", "EIGENVECTOR"), prob=0.1) {
  # Preserve the first graph
  original <- graph
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  # Save the start time
  start <- as.numeric(Sys.time())

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }
  seed <- c()
  # Until the criteria is fulfilled
  while (criteria_value < threshold) {
    # Use approximation algorithm for large graphs
    if (vcount(graph) > 10000) {
      centrality <- get_centrality_scores(graph, centrality_method=centrality_method, estimate=TRUE, cutoff=3)
    } else if (vcount(graph) > 5000) {
      centrality <- get_centrality_scores(graph, centrality_method=centrality_method, estimate=TRUE, cutoff=5)
    } else {
      centrality <- get_centrality_scores(graph, centrality_method=centrality_method)
    }
    df <- data.frame(name=V(graph)$name, centrality=centrality)
    influential <- df$name[df$centrality == max(df$centrality)][1]
    seed <- c(seed, influential)
    graph <- largest_component(delete.vertices(graph, V(graph)[V(graph)$name == influential]))
    if (vcount(graph) <= 2) {
      break
    }
    # Recalculate the new value for criteria to supply to loop condition
    if (criteria == "TIME") {
      # In case of TIME, the loop should end if the time threshold is up
      criteria_value <- Sys.time()
    } else if (criteria == "BUDGET") {
      # In case of BUDGET, the loop should end if the required budget threshold is met
      criteria_value <- length(seed)
  }
  end <- as.numeric (Sys.time())
  output <- NULL
  output$influential_nodes <- seed
  output$time <- (end - start)
  output$influence <- get_influence(original, V(original)[seed], test_method=test_method, lt_threshold=prob)
  output
}


#' Influence mining framework implementation based on Collective Influence method
#'
#' @name ci_influence_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @param prob probability threshold for INFLUENCE_LT. Default is 0.1
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
ci_influence_framework <- function(graph, criteria=c("TIME","BUDGET"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), prob=0.1) {
  # Preserve the first graph
  original <- graph
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  # Save the start time
  start <- as.numeric(Sys.time())

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }
  seed <- c()
  if (criteria == "BUDGET") {
    # If the criteria is budget, then calculate for all nodes and return top k nodes
    output <- collective_influence_influential(graph, budget=threshold, test_method=test_method)
    output$influential_nodes <- as.numeric(output$influential_nodes)
  } else {
    while (criteria_value < threshold) {
      ci <- sapply(V(graph), function(x) { collective_influence(graph, neighborhood_distance=2, x) })
      df <- data.frame(name=V(graph)$name, ci=ci)
      influential <- df$name[df$ci == max(df$ci)][1]
      seed <- c(seed, influential)
      graph <- largest_component(delete.vertices(graph, V(graph)[V(graph)$name == influential]))
      if (vcount(graph) <= 2) {
        break
      }
      if (criteria == "TIME") {
        # If the criteria is time, then keep picking nodes and calculate CI for highest degre node till the time runs out
        criteria_value <- Sys.time()
    }
    end <- as.numeric (Sys.time())
    output <- NULL
    output$influential_nodes <- seed
    output$time <- (end - start)
    output$influence <- get_influence(original, V(original)[seed], test_method=test_method, lt_threshold=prob)
  }
  output
}


#' Influence mining framework implementation based on Pagerank method
#'
#' @name pagerank_influence_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @param prob probability threshold for INFLUENCE_LT. Default is 0.1
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
pagerank_influence_framework <- function(graph, criteria=c("TIME","BUDGET"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), prob=0.1) {
  # Preserve the first graph
  original <- graph
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  # Save the start time
  start <- as.numeric(Sys.time())

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }
  seed <- c()
  if (criteria == "BUDGET") {
    # If the criteria is budget, then calculate for all nodes and return top k nodes
    output <- pagerank_influential(graph, budget=threshold, test_method=test_method)
    output$influential_nodes <- as.numeric(output$influential_nodes)
  } else {
    while (criteria_value < threshold) {
      pr <- page_rank(graph)$vector
      df <- data.frame(name=V(graph)$name, pr=pr)
      influential <- df$name[df$pr == max(df$pr)][1]
      seed <- c(seed, influential)
      graph <- largest_component(delete.vertices(graph, V(graph)[V(graph)$name == influential]))
      if (vcount(graph) <= 2) {
        break
      }
      if (criteria == "TIME") {
        # If the criteria is time, then keep picking nodes and calculate CI for highest degre node till the time runs out
        criteria_value <- Sys.time()
    }
    end <- as.numeric (Sys.time())
    output <- NULL
    output$influential_nodes <- seed
    output$time <- (end - start)
    output$influence <- get_influence(original, V(original)[seed], test_method=test_method, lt_threshold=prob)
  }
  output
}


#' Influence mining framework implementation based on K-core decomposition method
#'
#' @name coreness_influence_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @param prob probability threshold for INFLUENCE_LT. Default is 0.1
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
coreness_influence_framework <- function(graph, criteria=c("TIME","BUDGET"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), prob=0.1) {
  # Preserve the first graph
  original <- graph
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  # Save the start time
  start <- as.numeric(Sys.time())

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }
  seed <- c()
  if (criteria == "BUDGET") {
    # If the criteria is budget, then calculate for all nodes and return top k nodes
    output <- coreness_influential(graph, budget=threshold, test_method=test_method)
    output$influential_nodes <- as.numeric(output$influential_nodes)
  } else {
    while (criteria_value < threshold) {
      coreness <- graph.coreness(graph, mode="all")
      df <- data.frame(name=V(graph)$name, coreness=coreness)
      influential <- df$name[df$coreness == max(df$coreness)][1]
      seed <- c(seed, influential)
      graph <- largest_component(delete.vertices(graph, V(graph)[V(graph)$name == influential]))
      if (vcount(graph) <= 2) {
        break
      }
      if (criteria == "TIME") {
        # If the criteria is time, then keep picking nodes and calculate CI for highest degre node till the time runs out
        criteria_value <- Sys.time()
    }
    end <- as.numeric (Sys.time())
    output <- NULL
    output$influential_nodes <- seed
    output$time <- (end - start)
    output$influence <- get_influence(original, V(original)[seed], test_method=test_method, lt_threshold=prob)
  }
  output
}


#' Influence mining framework implementation based on CELF greedy method
#'
#' @name celf_greedy_influence_framework
#' @param graph is the igraph object
#' @param criteria is the stopping criteria for influence framework. Supported values are TIME, BUDGET and SPREAD
#' @param threshold is the limit value for criteria function. The execution steps as soon as criteria meets the threshold value. If the criteria is TIME, then threshold should be defined in seconds
#' @param test_method specifies the method to measure influence. Value "RESILIENCE" (number of total nodes REMOVED (NOT THE REMAINING ones as in original resilience function) from the graph); "INFLUENCE_IC" (see simulate_ic method); "INFLUENCE_LT" (see simulate_lt method). Default is "RESILIENCE"
#' @param prob probability threshold for INFLUENCE_LT. Default is 0.1
#' @return output object containing influential nodes, time taken and the influence of nodes according the test_method
#' @import igraph
#' @export
celf_greedy_influence_framework <- function(graph, criteria=c("TIME","BUDGET"), threshold, test_method=c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC"), prob=0.1) {
  seed <- c()
  # For each node, compute its influence independently and add maintain in a table, while setting the flag=0
  V(graph)$name <- 1:vcount(graph)
  df <- data.frame(node=V(graph)$name, gain=0, flag=0)
  # Save the start time
  start <- as.numeric(Sys.time())

  for (i in 1:vcount(graph)) {
    df$gain[i] <- get_influence(graph, V(graph)[i], test_method=test_method, lt_threshold=prob)
  }
  # Arrange the data frame by marginal gains
  df <- arrange(df, desc(gain))

  # Modification for the framework
  if (criteria == "TIME") {
    # If the criteria is TIME, then the threshold should be set according to time
    criteria_value <- Sys.time()
    threshold <- criteria_value + seconds(threshold)
  } else {
    criteria_value <- 0
  }

  # Until the criteria is fulfilled
  while (criteria_value < threshold) {
    top_row <- df[1,]
    u <- V(graph)[top_row$node]
    # If the flag is the size of the seed set so far, then add this node to the seed and remove from data frame
    if (top_row$flag == length(seed)) {
      seed <- c(seed, top_row$node)
      df <- df[-1,]
      if (nrow(df) == 0) {
        break
      }
    } else {
      # Otherwise compute the marginal gain with this node
      current_influence <- get_influence(graph, V(graph)[seed], test_method=test_method, lt_threshold=prob)
      top_row$gain <- get_influence(graph, V(graph)[c(seed, u)], test_method=test_method, lt_threshold=prob) - current_influence
      # Store the length of seed in the flag
      top_row$flag <- length(seed)
      # Update the values for this row in data frame
      df[1,] <- top_row
      # Sort the data frame again by gain
      df <- arrange(df, desc(gain))
    }
    # Recalculate the new value for criteria to supply to loop condition
    if (criteria == "TIME") {
      # In case of TIME, the loop should end if the time threshold is up
      criteria_value <- Sys.time()
    } else if (criteria == "BUDGET") {
      # In case of BUDGET, the loop should end if the required budget threshold is met
      criteria_value <- length(seed)
  }
  end <- as.numeric (Sys.time())
  output <- NULL
  output$influential_nodes <- V(graph)[seed]$name
  output$time <- (end - start)
  output$influence <- get_influence(graph, output$influential_nodes, test_method=test_method, lt_threshold=prob)
  output
}



#####################################
##### COMPREHENSIVE BUDGET TEST #####
#####################################
graphs <- load_graphs()

# TEST adaptive centrality methods against BUDGET criteria
output_file <- paste0(root_dir, "influence_framework_budget.csv")
criteria <- "BUDGET"
thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05)
test_methods <- c("RESILIENCE", "INFLUENCE_LT", "INFLUENCE_IC")
graphs <- list(karate, football, world_trade, nematode, automata, author, protein, influence_citation, routers, hamster, moreno, fb_tvshows, political, usairport, city2010, city2013, city2016, city2019, spam, arxiv, ita2000, jdk6)

plan(multiprocess, workers=4)
future_lapply(graphs, function(graph) {
  for (test_method in test_methods) {
    for (t in thresholds) {
      threshold <- ceiling(vcount(graph) * t)
      # Degree
      centrality_method <- "DEGREE"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_degree", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Betweenness
      centrality_method <- "BETWEENNESS"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_betweenness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Closeness
      centrality_method <- "CLOSENESS"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_closeness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Eigenvector
      centrality_method <- "EIGENVECTOR"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_eigenvector", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Pagerank
      output <- pagerank_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="pagerank", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Coreness
      output <- coreness_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="coreness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Collective influence
      output <- ci_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="collective_influence", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # CELF Greedy
      output <- celf_greedy_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="celf_greedy", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')
    }
  }
})

# Large graphs
thresholds <- c(0.00125)
graphs <- list(johnshopkins, caida, linux, youtube)

plan(multiprocess, workers=4)
future_lapply(graphs, function(graph) {
  for (test_method in test_methods) {
    for (t in thresholds) {
      threshold <- ceiling(vcount(graph) * t)
      # Degree
      centrality_method <- "DEGREE"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_degree", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Closeness
      centrality_method <- "CLOSENESS"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_closeness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Eigenvector
      centrality_method <- "EIGENVECTOR"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_eigenvector", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Pagerank
      output <- pagerank_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="pagerank", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Coreness
      output <- coreness_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="coreness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Betweenness
      centrality_method <- "BETWEENNESS"
      output <- adaptive_centrality_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method, centrality_method=centrality_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="adaptive_betweenness", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # Collective influence
      output <- ci_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="collective_influence", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')

      # CELF Greedy
      output <- celf_greedy_influence_framework(graph=graph, criteria=criteria, threshold=threshold, test_method=test_method)
      results <- data.frame(graph=graph$name, test_method=test_method, criteria=criteria, size=vcount(graph), threshold=t, influence_function="celf_greedy", time=output$time, influence=output$influence, budget=length(output$influential_nodes), influential_nodes=toString(toJSON(unlist(output$influential_nodes))))
      print(results)
      write.table(results, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')
    }
  }
})


# Load in Workbook
file <- "results/consolidated_results.csv"
results <- read.csv2(file, header=TRUE, sep=',', stringsAsFactors=FALSE)
results$threshold <- as.numeric(results$threshold)
results$time <- round(as.numeric(results$time), 6)

budget <- results[results$criteria == "BUDGET", ]
networks <- unique(budget$graph)
methods <- unique(budget$influence_function)


###########################################
### Analyze results from BUDGET results ###
###########################################

criteria <- "BUDGET"
x <- results[results$criteria == criteria, ]
summarized <- data.frame(criteria=c(), test_method=c(), threshold=c(), graph=c(), max_inf=c(), max_inf_function=c())
thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05)

for (test_method in c("RESILIENCE", "INFLUENCE_IC", "INFLUENCE_LT")) {
  # In each graph
  for (graph in unique(x$graph)) {
    # And threshold
    for (threshold in thresholds) {
      # Print the algorithm with highest influence
      y <- x[x$criteria == criteria & x$test_method == test_method & x$threshold == threshold & x$graph == graph,]
      if (nrow(y) == 0) {
        threshold <- 0.00125
        y <- x[x$criteria == criteria & x$test_method == test_method & x$threshold == threshold & x$graph == graph,]
      }
      max_inf <- max(y$influence)
      max_inf_function <- y$influence_function[y$influence == max_inf]
      print(paste(criteria, test_method, threshold, graph))
      df <- data.frame(criteria=c(), test_method=c(), threshold=c(), graph=c(), max_inf=c(), max_inf_function=c())
      for (mif in max_inf_function) {
        df <- rbind(df, data.frame(criteria, test_method, threshold, graph, max_inf, mif))
      }
      print(graph)
      write.table(df, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')
      if (!is.na(max_inf))
        summarized <- rbind(summarized, df)
    }
  }
}

summarized
# Copy to clipboard to move to Excel
clipr::write_clip(summarized)


x <- results[results$criteria == criteria & results$graph %in% c('fb_johnshopkins', 'caida', 'linux_code'), ]
summarized <- data.frame(criteria=c(), test_method=c(), threshold=c(), graph=c(), min_time=c(), min_time_function=c())
thresholds <- c(0.005, 0.01, 0.025)

for (test_method in c("RESILIENCE", "INFLUENCE_IC", "INFLUENCE_LT")) {
  # In each graph
  for (graph in unique(x$graph)) {
    # And threshold
    for (threshold in thresholds) {
      print(paste(graph, test_method))
      print(threshold)
      # Print the algorithm with highest influence
      y <- x[x$criteria == criteria & x$test_method == test_method & x$threshold == threshold & x$graph == graph,]
      if(nrow(y) == 0) {
        print("FAILED! FAILED! FAILED!")
        next
      }
      min_time <- min(y$time)
      min_time_function <- y$influence_function[y$time == min_time]
      print(paste(criteria, test_method, threshold, graph))
      df <- data.frame(criteria=c(), test_method=c(), threshold=c(), graph=c(), min_time=c(), min_time_function=c())
      for (mtf in min_time_function) {
        df <- rbind(df, data.frame(criteria, test_method, threshold, graph, min_time, mtf))
      }
      write.table(df, output_file, row.names=FALSE, col.names=FALSE, append=TRUE, sep=',')
      if (!is.na(max_inf))
        summarized <- rbind(summarized, df)
    }
  }
}

summarized
# Copy to clipboard to move to Excel
clipr::write_clip(summarized)
