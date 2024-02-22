---
title: No title

institute:
  - stat: Department of Statistics, Lund University
  - medfak: Division of Infection Medicine, Faculty of Medicine, Lund University
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
  - \usepackage{algorithm2e}

---

# Some text

Hello

# Algorithm 1

Just a sample algorithm
\begin{algorithm}[H]
\DontPrintSemicolon
\SetAlgoLined
\KwResult{Write here the result}
\SetKwInOut{Input}{Input}\SetKwInOut{Output}{Output}
\Input{Write here the input}
\Output{Write here the output}
\BlankLine
\While{While condition}{
    instructions\;
    \eIf{condition}{
        instructions1\;
        instructions2\;
    }{
        instructions3\;
    }
}
\caption{While loop with If/Else condition}
\end{algorithm} 


