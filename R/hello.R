# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {print("Hello, world!")}

# 가중치 구하는 함수
get.weigt <- function(x) {
  # 순위 자체가 가중치
  rank <- function(x) {
    x
  }

  # 서열합 가중치(rank-sum weight)
  rank.sum <- function(x) {
    val <- max(x) - x + 1
    val/sum(val)
  }
  # 서열 역수 가중치(rank reciprocal weight)
  rank.reci <- function(x) {
    inver <- 1/x
    inver/sum(inver)
  }
  # 서열 중심법 ROC(rank order centrord)
  roc <- function(x) {
    sort.x <- sort(x)
    n <- length(x)
    sapply(x, function(x) sum(1/sort.x[x:n])/n)
  }

  methods <- list(rank = rank,
                  rank.sum = rank.sum,
                  rank.reci = rank.reci,
                  roc = roc)
  methods[[x]]
} ## 가중치 계산 함수(end)

# weighted
## 정규화 함수
normalize <- function(x) {
  # max
  linear.max <- function(x) {
    apply(x, 2, function(x) x/max(x))
  }
  # min
  linear.sum <- function(x) {
    apply(x, 2, function(x) x/sum(x))
  }
  # mean
  linear.mean <- function(x) {
    apply(x, 2, function(x) x/mean(x))
  }
  # proportional score
  prop <- function(x) {
    apply(x, 2, function(x) 1-(max(x)-x)/(max(x)-min(x)))
  }
  # 벡터 정규화
  norm.vec <- function(x) {
    square.sums <- apply(x^2, 2, sum)
    t(t(x)/sqrt(square.sums))
  }

  methods <- list(linear.max = linear.max,
                  linear.sum = linear.sum,
                  linear.mean = linear.mean,
                  prop = prop,
                  norm.vec = norm.vec)
  methods[[x]]
} ## 정규화 함수(end)

## 주고유 벡터 구하기
get.Peigen <- function(x) {
  prin.eigen <- function(x) {
    e <- eigen(x)
    p <- e$vector[,1]
    p
  }

  row.sums <- function(x) {
    ri <- apply(x, 1, sum)
    a <- ri / sum(ri)
    a
  }

  col.sums <- function(x) {
    cj <- 1 / (apply(x, 2, sum))
    b <- cj / sum(cj)
    b
  }

  Brow.sums <- function(x) {
    cj <- apply(x, 2, sum)
    bi <- t(apply(x, 1, function(x) x/cj))
    b.sums <- apply(bi, 1, sum)
    b.sums/length(cj)
  }

  geo.mean <- function(x) {
    p  <- apply(x,1, function(x) Reduce("*",x))
    pi <- p ^ (1/length(p))
    x  <- pi / sum(pi)
    x
  }
  methods <- list(prin.eigen = prin.eigen,
                  row.sums = row.sums,
                  col.sums = col.sums,
                  Brow.sums = Brow.sums,
                  geo.mean = geo.mean)
  methods[[x]]
} ##고유벡터 함수(end)


## 효율전선 그리기
#
sign <-function(p1, p2, p3) {
  p1.x <- p1[1]
  p1.y <- p1[2]
  p2.x <- p2[1]
  p2.y <- p2[2]
  p3.x <- p3[1]
  p3.y <- p3[2]

  return ((p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y))
}

# 세 점 안의 점 찾기
PointInTriangle <- function(pt, v1, v2, v3) {
  pt <- pt
  v1 <- v1
  v2 <- v2
  v3 <- v3

  d1 = sign(pt, v1, v2);
  d2 = sign(pt, v2, v3);
  d3 = sign(pt, v3, v1);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  dodo <- all(pt == v1) | all(pt == v2) | all(pt == v3)
  if( dodo ) res <- F else res <- !(has_neg && has_pos)
  return (res)
}

# 두 점의 거리구하기
distance <- function(a,b) {
  a.x <- a[1]
  a.y <- a[2]
  b.x <- b[1]
  b.y <- b[2]
  return(sqrt((a.x - b.x)**2 + (a.y - b.y)**2))
}

# 직선위의 점 찾기
PointOnLine <- function(pt, v1, v2) {
  return(distance(v1,pt) + distance(pt,v2) == distance(v1,v2))
}

## 입력값 설명
#1. alts : 대안 열(열 번호 또는 열 이름(문자형 X))
#2. data : 대안을 포함은 원본데이터
#3. rank : 평가기준 선호 순위(정수)
#4. cost : 비용함수 골라냄 T,'-' : cost 함수
#5. eff.cost.ratio : 효익 비용 가중치 직접 설정
#6 linear : 비례 선형번환(T) default : F

smart <- function(alts, data, rank, cost, get.w = "roc", eff.cost.ratio = c(1, 1), linear = F) {
  alts_name <- substitute(alts)
  data_name <- substitute(data)
  rank_name <- substitute(rank)
  cost_name <- substitute(cost)
  eff.cost.ratio_name <- substitute(eff.cost.ratio)
  linear_name <- substitute(linear)

  # alts를 열 번호로 변환
  ## error 코드
  alts.char <- deparse(alts_name)
  if ( length(alts.char) != 1 ) {
    return(simpleError("argument alts should be length one vector", alts_name))
  } else if (!any(names(data) == alts.char)) {
    return(simpleError("There is no such column in given data", alts_name))
  }
  # alts.num : 열 번호
  nl <- as.list(seq_along(data))
  names(nl) <- names(data)
  alts.num <- eval(substitute(alts), nl, parent.frame())
  #alts.num <- 5
  # 대안 변수 : 제거 및 행 이름으로 변환
  alts.name <- as.character(data[,alts.num]) # 대안 열 이름
  dd <- data[-alts.num]
  row.names(dd) <- alts.name
  criteria.name <- colnames(dd) # 평가기준 열 이름

  orin.dd <- dd # 기존 data 저장

  ## cost 변수 역수치
  # cost 오류 및 논리형 변환
  if( !missing(cost) ) {
    if (ncol(dd) != length(cost)) {
      return(simpleError("wrong length argument to cost", cost_name))
    } else if ( !is.character(cost) & !is.logical(cost) ) {
      return(simpleError("aurgument cost should be a logical or character vector('+' or '-')"))
    } else if( is.character(cost) ){
      if ( !all( cost == '+' | cost == '-') ) {
        return(simpleError("aurgument cost should be a logical or character vector('+' or '-')", cost_name))
      } else cost <- ifelse(cost == '+', F, T)
    }
  } else cost <- rep(F, ncol(dd))


  ## 비례 선형번환
  if( linear == T ) {
    linear.tran <- function(x, y) {
      val <- (max(x)- x)/(max(x)-min(x))*100
      if( y == T ) val else 100 - val
    }
    dd <- as.data.frame(Map(linear.tran, dd, cost))
    row.names(dd) <- alts.name
    linear.dd <- dd # 비례 선형번환 값 저장
  } else linear.dd <- NULL

  if(max(dd) > 100 | min(dd) < 0 ) return(simpleError("elements should be in range max 100, min 0 OR give TRUE to argument linear"))

  w.fun <- get.weigt(get.w) # w.fun 지정된 가중치 함수

  # 효익
  eff.col <- colnames(dd)[!cost]
  eff.val <- dd[!cost]
  eff.w <- w.fun(order(rank[!cost]))

  # 비용
  cost.col <- colnames(dd)[cost]
  cost.val <- dd[cost]
  cost.w <- w.fun(order(rank[cost]))

  # 종합 효익 가치
  eff.score <- as.matrix(eff.val)%*%eff.w
  rownames(eff.score) <- alts.name
  cost.score <- as.matrix(cost.val)%*%cost.w
  rownames(cost.score) <- alts.name

  # 대안 점수
  eff.cost <- eff.cost.ratio/sum(eff.cost.ratio)
  scores <- eff.score*eff.cost[1] + cost.score*eff.cost[2]
  rownames(scores) <- alts.name

  res <- list(data = orin.dd, # 원본데이터
              criteria.name = criteria.name, # 평가 기준 열 이름
              alts.col_name = alts_name, # 대안 열 이름
              alts.name = alts.name, # 대안 항목명
              linear = linear, # 비례 선형번환 여부
              linear.data = linear.dd, # 비례 선형번환 data
              get.w = get.w, # ?

              cost = cost, # 비용 함수 여부

              eff.col = eff.col, # 효익 열 이름
              eff.w = eff.w, # 효익 가중치
              eff.score = eff.score, # 효익 점수

              cost.col = cost.col, # 비용 열 이름
              cost.w = cost.w, # 비용 가중치
              cost.score = cost.score, # 비용 점수

              eff.cost = eff.cost, # 효익 비용 점수
              scores = scores # 대안 점수
  )
  class(res) <- "smart"
  return(res)
}
# 기본형
res <- smart(alts = company, data = dat, rank = c(1,3,2,4), cost = c(F,T,F,F), linear = F)
res

# 비례점수 적용
res <- smart(alts = company, data = dat2, rank = c(1,3,2,4), cost = c(F,T,F,F), linear = T, eff.cost.ratio = c(30, 70))
res

# 정규화 되지 않은 data
# res <- smart(alts = company, data = dat2, rank = c(1,3,2,4), cost = c(F,T,F,F), linear = F, eff.cost.ratio = c(30, 70))
# res  # error 발생

summary.smart <- function(x) {
  out.mat <- NULL
  cat('\t\tSMART\n\n')
  cat('Alternative column :', x$alts.col_name, '\n')
  cat('\n')

  weight <- x$eff.cost.ratio
  if( x$linear ) {
    cat("*proportional score\n")
    data <- x$linear.data
  } else data <- x$data

  if( any(!x$cost) ) {
    cat("< Efficiency :", x$eff.cost[2], ">\n")
    eff.mat <- rbind(data[!x$cost], x$eff.w)
    eff.mat <- cbind(eff.mat, c(x$eff.score, NA))
    row.names(eff.mat)[nrow(eff.mat)] <- "eff.Weight"
    colnames(eff.mat)[ncol(eff.mat)] <- "eff.Score"
    eff.mat <- round(data.frame(eff.mat), 2)
    eff.mat[nrow(eff.mat),ncol(eff.mat)] <- "-"
    print(eff.mat)
    cat("\n")
  }
  if( any(x$cost) ) {
    cat("< Cost :", x$eff.cost[1], ">\n")
    cost.mat <- rbind(data[x$cost], x$cost.w)
    cost.mat <- cbind(cost.mat, c(x$cost.score, NA))
    row.names(cost.mat)[nrow(cost.mat)] <- "cost.Weight"
    colnames(cost.mat)[ncol(cost.mat)] <- "cost.Score"
    cost.mat <- round(data.frame(cost.mat), 2)
    cost.mat[nrow(cost.mat),ncol(cost.mat)] <- "-"
    print(cost.mat)
    cat("\n")
  }
  cat('method :', x$get.w, '\n')
  cat('\n')


  scores <- t(data.frame(x$scores))
  colnames(scores) <- x$alts.name
  row.names(scores) <- "Scores"
  print(round(scores,4)); cat('\n')

  sorted.score <- sort(-scores)
  s.names <- x$alts.name[order(-scores)]
  ranking <- s.names[1]

  for (i in 2:length(scores)) {
    equal <- sorted.score[i-1] == sorted.score[i]
    ranking <- if ( equal ) {
      paste(ranking, "=", s.names[i])
    } else paste(ranking, "", s.names[i])
  }

  cat('Alternatives Ranking :', ranking)
}
summary.smart(res)

