# Math::Regressionizer

Raku package with a class that allows pipeline-like specification and execution of regression workflows.

Extensive guide is given in the Jupyter notebook 
["Rapid-specification-of-regression-workflows.ipynb"](https://github.com/antononcube/Raku-Math-Regressionizer/blob/main/docs/Rapid-specification-of-regression-workflows.ipynb)
and the corresponding Markdown document
["Rapid-specification-of-regression-workflows.md"](https://github.com/antononcube/Raku-Math-Regressionizer/blob/main/docs/Rapid-specification-of-regression-workflows.md).

------

## Features summary 

- The class `Math::Regressionizer` facilitates rapid specifications of regressions workflows.
  - To quickly specify: 
    - data rescaling and summary
    - regression computations
    - outliers finding
    - conditional Cumulative Distribution Functions (CDFs) reconstruction
    - plotting of data, fits, residual errors, outliers, CDFs 

- `Math::Regressionizer` works with datasets, numerical arrays, lists of numbers, and lists of numeric pairs.


### Details and arguments

- The curves computed with Quantile Regression are called **regression quantiles**.

- `Math::Regressionizer` has three regression methods:
  - `quantile-regression`
  - `quantile-regression-fit`
  - `least-squares-fit`
  
- The regression quantiles computed with the methods `quantile-regression` and `quantile-regression_fit` 
  correspond to probabilities specified with the argument `:probs(:@probabilities)`.

- The method`quantile-regression` computes fits using a B-spline functions basis.  
  - The basis is specified with the arguments `:@knots` and `:order(:$interpolation-order)`.
  - `:order(:$interpolation-order)` is 3 by default. 

- The methods `quantile-regession-fit` and `least-squares-fit` use lists of basis functions to fit with 
  specified with the argument `:funcs(:@basis)`.

### Workflows flowchart

The following flowchart summarizes the workflows that are supported by `Regressionizer`:

![](https://raw.githubusercontent.com/antononcube/Python-Regressionizer/main/docs/img/Quantile-regression-workflow-extended.jpg)


------

## Usage examples 

Load packages:

```raku
use Math::Regressionizer;
use Data::Generators;
```

Generate random data:

```raku
srand(0);
my @x = 0, 0 + 1/300 ... 2;
my @y = @x.map({ sin(2 * pi * $_) }) <<+>> random-variate(NormalDistribution.new(0, 0.4), @x.elems);
my @data = @x Z @y;
```

Compute quantile regression for probabilities `[0.2, 0.5, 0.8]` and make the corresponding plot:

```raku
#% html
my $obj = 
  Math::Regressionizer.new(@data)
  .quantile-regression(knots => 8, probs => [0.2, 0.5, 0.8])
  .plot(title => "B-splines fit", template => "dark", format => 'html')
```

**Remark** If the code cell above is executed in a Jupyter notebook with a Raku kernel then the plot is shown.
"Math::Regressionizer" uses the package ["JavaScript::Google::Charts"](https://github.com/antononcube/Raku-JavaScript-Google-Charts) for its plots.

Export the plot obtained above:

```raku
spurt('qr-plot.html', $obj.take-value)
```

![](https://raw.githubusercontent.com/antononcube/Python-Regressionizer/main/docs/img/random-data-B-spline-rqs.png)

------

## References

### Articles, books

[RK1] Roger Koenker, 
[Quantile Regression](https://books.google.com/books/about/Quantile_Regression.html?id=hdkt7V4NXsgC), 
Cambridge University Press, 2005.

[RK2] Roger Koenker,
["Quantile Regression in R: a vignette"](https://cran.r-project.org/web/packages/quantreg/vignettes/rq.pdf),
(2006),
[CRAN](https://cran.r-project.org/).

[AA1] Anton Antonov,
["A monad for Quantile Regression workflows"](https://github.com/antononcube/MathematicaForPrediction/blob/master/MarkdownDocuments/A-monad-for-Quantile-Regression-workflows.md),
(2018),
[MathematicaForPrediction at GitHub](https://github.com/antononcube/MathematicaForPrediction).

### Packages, paclets

[RKp1] Roger Koenker,
[`quantreg`](https://cran.r-project.org/web/packages/quantreg/index.html),
[CRAN](https://cran.r-project.org/).

[AAp1] Anton Antonov,
[Quantile Regression, WL paclet](https://github.com/antononcube/WL-QuantileRegression-paclet),
(2014-2023),
[GitHub/antononcube](https://github.com/antononcube).

[AAp2] Anton Antonov,
[Monadic Quantile Regression, WL paclet](https://github.com/antononcube/WL-MonadicQuantileRegression-paclet),
(2018-2024),
[GitHub/antononcube](https://github.com/antononcube).

[AAp3] Anton Antonov,
[`QuantileRegression`](https://resources.wolframcloud.com/FunctionRepository/resources/QuantileRegression),
(2019),
[Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/resources/QuantileRegression).

### Repositories

[AAr1] Anton Antonov,
[DSL::English::QuantileRegressionWorkflows, Raku package](https://github.com/antononcube/Raku-DSL-English-QuantileRegressionWorkflows),
(2020),
[GitHub/antononcube](https://github.com/antononcube/Raku-DSL-English-QuantileRegressionWorkflows).

### Videos

[AAv1] Anton Antonov,
["Boston useR! QuantileRegression Workflows 2019-04-18"](https://www.youtube.com/watch?v=a_Dk25xarvE),
(2019),
[Anton Antonov at YouTube](https://www.youtube.com/@AAA4Prediction).

[AAv2] Anton Antonov,
["useR! 2020: How to simplify Machine Learning workflows specifications"](https://www.youtube.com/watch?v=b9Uu7gRF5KY),
(2020),
[R Consortium at YouTube](https://www.youtube.com/channel/UC_R5smHVXRYGhZYDJsnXTwg).

[AAv3] Anton Antonov,
["Quantile Regression Workflows (PyData Global 2024)"](https://www.youtube.com/watch?v=mYXCynTPHsA),
(2025),
[Anton Antonov at YouTube](https://www.youtube.com/@AAA4Prediction).
