---
title: No title

institute:
  - medfak: Division of Infection Medicine, Faculty of Medicine, Lund University
  - stat: Department of Statistics, Lund University

author: 
  - Erik Hartman:
      institute: medfak
  - Jonas Wallin:
      institute: stat

output:
  pdf_document:
    latex_engine: xelatex
    pandoc_args: ["--lua-filter=scholarly-metadata.lua","--lua-filter=author-info-blocks.lua"]

urlcolor: "blue"
bibliography: bibliography.bib

header-includes:
  - \usepackage[ruled,vlined,linesnumbered]{algorithm2e}

---


## Gradient

$V$ är alla peptider.
$v$ är en peptid.

$ \Gamma = (V,E) $ är den riktade grafen.
$ P (\cdot | v) $ är en fördelning av peptider (vektor av storlek $m$ som är antalet unika peptider) givet att du står i peptid $v$.

$V_L$ är den längsta peptiden (nod utan in-kanter).

$P(v_L \rightarrow v)$ är sannolikheten att besöka $v$. Också samma som $P(v_L \rightarrow v) = P(v | v_L) / w_v^v$

$ P(\cdot | v) = \bold1_v * w^v_v + \sum_{i \sim v}P(\cdot |i)*w^v_i$

$ P(AB) = \bold1_{AB} * w_{AB}^{AB} + P(A)*w_A^{AB} + P(B)*w_B^{AB}$

$w_v^v = P(\text{stanna i v})$ därför är sannolikhetsfördelningen av peptider, $P(\cdot)$ en vektor.

Notera att $w_v^v = 1 - \sum_{i \sim v}{w^i_v}$

Därför blir gradienten:
$ \frac{dP(\cdot | v_L)}{dw^i_v} = P(v_L \rightarrow v)(P(\cdot|i)-\bold1_v)$

$D_{KL} = log(P(\cdot|v_L)) = \sum_{v \in V}log(P(v|v_L))$

$n_v$ är antalet observerade peptider av typ $v$. 

$\frac{dD_{KL}}{dw^v_i} = \frac{d}{dw^v_i}(P(\cdot | v_L))^T * ( \frac{n_i}{P(v_i | v_L)} + ... \frac{n_m}{P(v_m | v_L)})$ där indexering av peptider är i ordningen av vektorn $P(\cdot | v_L)$. Alternativt kan skrivas med summa.

Constraints: 



## Algorithms


\begin{algorithm}[H]
\caption{Simulate proteolysis of protein $P$.}
\SetAlgoLined
\DontPrintSemicolon
\KwIn{protein $P$, $n_{generate}$ $\theta_{enzyme}$, $\theta_{gamma}$, $p_{endo}$, $p_{exo}$}
\KwOut{$f$}

$f \gets $ dict[sequence : count]\;
$f(P) \gets n_P$\;

\While{$i < n_{generate}$}{
    $n_{generate} \gets \sum_xT(x)$\;
    $x \sim U(0,1) $\;
    \eIf{$x>p_{endo}$}{
        //exoprotease\;
        sequence to chew $s \sim f$ //sample sequence from dict with weights count\;
        $x \sim U(0,1) $\;
        $a \gets gamma(len(s), \theta_{gamma})$\;
        \If{$a<x$}{
            //accept\;
            $F(P) \gets F(P)+1$\;
            $n_{generate} \gets n_{generate} +1$\;
        }
    }
    {
        //endoprotease\;
        $f_{cut}(s) \gets \sum_s N_{aa}^{s}*\theta_{aa}$\;
        sequence to cut $s \gets f_{cut}(s)$//sequence to cut \;
        first index to cut $index_a \sim s(\theta_{aa}(aa_x))$\;
        $index_2 \sim s(\theta_{aa}(aa_x))*gamma(|{index_1 - index_2}|)$\;
        left $\gets s[: min(index_1, index_2)+1]$\;
        middle $\gets s[min(index_1, index_2)+1:max(index_1, index_2)]$\;
        right $\gets s[max(index_1, index_2)+1:]$\;
        \If{len(middle)>5}{
            $f(middle) \gets f(middle) + 1$\;
            $n_{generate} \gets n_{generate} +1$\;
        }
        \For{$s$ in [left, right]}{
            $x \sim U(0,1) $\;
            $a \gets gamma(len(s), \theta_{gamma})$\;
            \If{$a<x$ and $s$ is not terminal peptide in $P$}{
                //accept\;
                $F(P) \gets F(P)+1$\;
                $n_{generate} \gets n_{generate} +1$\;
            }
        }
    }
}

Return $f$
\end{algorithm}


\begin{algorithm}[H]
\caption{Estimating $\theta$ numerically. To generate a guess, simulate degradation of protein $P$ with parameters $\theta$ to generate $n_{generate}$ peptides (see Algorithm 1).}
\SetAlgoLined
\KwIn{protein $P$, $n_P$, true distribution $T$, $\theta$, $lr_{endo}$, $lr_{exo}$}
\KwOut{$\theta$}


\For{i from 0 to $n_{endo}$}{
    Generate guess $G$\;
    Compute loss $L \gets D_{KL}(G||T) + D_{KL}(T||G)$\;
    \For{each amino acid aa in $\theta_{aa}$}{
        $\theta_{aa}(aa) \gets \theta_{aa}(aa) + lr_{endo}$\;
        Generate guess $\hat{G}$ with new $\theta$\;
        Compute new loss $\hat{L} \gets D_{KL}(\hat{G}||T) + D_{KL}(T||\hat{G})$\;
        \While{$\hat{L} < L$}{
            Compute weighted learning rate $ lr_w \gets lr_{endo}\;*\hat{L}-L$\;
            $L \gets \hat{L}$\;
            $\theta_{aa}(aa) \gets \theta_{aa}(aa) + lr_w$\;
            Generate guess $\hat{G}$ with new $\theta$\;
            Compute new loss $\hat{L} \gets D_{KL}(\hat{G}||T) + D_{KL}(T||\hat{G})$\;
        }
       $\theta_{aa}(aa) \gets \theta_{aa}(aa) - lr_{endo}$ //revert the initial parameter-change (before while-loop)\;
    }
}

Generate guess $G$\;
Compute loss $L \gets D_{KL}(G||T) + D_{KL}(T||G)$\;

\For{i from 0 to $n_{exo}$}{
    $x \gets$ uniformly random from ${1, -1}$\;
    $e \gets lr_{exo} * x$\;
    $\theta_{exo} \gets \theta_{exo} + e$\;
    Generate guess $\hat{G}$ with new $\theta$\;
    Compute new loss $\hat{L} \gets D_{KL}(\hat{G}||T) + D_{KL}(T||\hat{G})$\;
    \eIf{$\hat{L} > L$}{
        $\theta_{exo} \gets \theta_{exo} - e$\;
    }{
        $L \gets \hat{L}$\;
    }
}
Return $\theta$\;

\end{algorithm}



\begin{algorithm}[H]
\caption{Estimate weights in graph with gradient descent.}
\SetAlgoLined
\DontPrintSemicolon
\KwIn{observed distribution $T$, $lr$, $\lambda_1$, $\lambda_2$, $n_{iterations}$}
\KwOut{$G$}

$G \gets \{V,E,W\}$\;
$V_L$ denotes the original protein (bottom node)\;



\For{i from 0 to $n_{iterations}$}{
    \;
    //generate guess\;
    $p_{generated} \gets \{\}$// this represents the output distribution if starting from a given node\;
    $T = \{ s \in V \mid \text{there exists no } (s, t) \in E \text{ for any } t \in V \}$ //terminal nodes\;
    \For{$t \in T$}{
        $p_{generated}[t] \gets \bold{1}_t$ //onehot\;
    }
    \While{$\text{all nodes in } G \text{ is not solved}$}{
        solvable $\gets \{s \in V \mid t \in p_{generated} \text{ for all } (s,t) \in E \}$\;
        \For{$s \in \text{solvable}$}{
            $p_{generated}[s] = \sum{w_{s,t}*p_{generated}[t]} + 1-\sum{w_{s,t}} *\bold1_t$\;
        }
    }
    $\hat{T} \gets p_{generated}[V_L]$\;
    \;
    //compute loss\;
    $L_1 \gets \lambda_1\sum{|w|}$\;
    $L_2 \gets \lambda_2\sum{w^2}$\;
    $L \gets D_{KL}(T \mid \hat{T}) + D_{KL}(\hat{T} \mid T) +L_1 + L_2$\;
    \;
    //compute gradient\;
    $\frac{dT}{dw} \gets \hat{T}_{V_L, t}(\hat{T}_s - \bold1_s)$\;
    $\frac{dL}{dT} \gets -\frac{T}{\hat{T}}$\;
    $\frac{dL}{dw} \gets \frac{dL}{dT}*\frac{dT}{dw}$\;
    \;
    //update graph\;
    $k \gets 1$\;
    \For{$s \in V$}{
        $\hat{w}_{s,t} \gets max(0, w_{s,t} - lr*(\frac{dL}{dw})_{s,t})$\;
        $d_{s,t} = \hat{w}_{s,t} - w_{s,t}$\;
    }
    
    \While{ $\sum_t{w_{s,t} + k*d_{s,t}} > 1$ }{
        $k \gets k / 2$ \;
    }

    \While{$\text{a better graph is not found or $k$ isn't extremely small}$}{
        $\hat{W} \gets W + d*k$\;
        $\hat{T} \gets \text{ generate guess with new weights}$\;
        $\hat{L} \gets D_{KL}(T \mid \hat{T}) + D_{KL}(\hat{T} \mid T) +L_1 + L_2$\;
        \eIf{$\hat{L} \le L$}{
            $G \gets \{V,E,\hat{W}\}$\;
        }{
            $k \gets k/2$\;
        }
    }
}
Return $G$
\end{algorithm}


