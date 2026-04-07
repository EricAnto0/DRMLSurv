# Internal helper to annotate errors with a step label

Internal helper to annotate errors with a step label

## Usage

``` r
tag_error(expr, label)
```

## Arguments

- expr:

  Expression to evaluate.

- label:

  Character scalar used to prefix any thrown error.

## Value

The result of `expr`, or an error with a labeled message.
