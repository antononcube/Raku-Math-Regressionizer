use v6.d;

unit class Math::Regressionizer;

use Math::QuantileRegression;
use JavaScript::Google::Charts;

my &core-quantile-regression = &quantile-regression;
my &core-quantile-regression-fit = &quantile-regression-fit;

# A small, stateful façade over Math::QuantileRegression.
# Methods return the invocant so that a regression workflow can be written as a pipeline.
has @!data;
has @!basis-functions;
has %!regression-quantiles;
has $!value;

#==========================================================
# Creators
#==========================================================

multi method new($data) {
    my $object = self.bless;
    $object.set-data($data);
    $object
}

#==========================================================
# Verifiers
#==========================================================

method !probabilities($probs --> Array:D) {
    my @probs = $probs.defined ?? ($probs ~~ Positional:D ?? $probs.List !! [$probs]) !! [0.25, 0.5, 0.75];
    die 'Probabilities must be numeric values between zero and one.'
    unless @probs.elems && @probs.all ~~ Numeric:D && @probs.all >= 0 && @probs.all <= 1;
    return @probs.Array
}

method !require-quantiles() {
    die 'Cannot find regression functions; compute a regression first.' unless %!regression-quantiles.elems;
}

#==========================================================
# Setters / takers
#==========================================================

method set-data($data --> ::?CLASS:D) {
    die 'Data must be a non-empty positional collection.'
    unless $data ~~ Positional:D && $data.elems;

    @!data = ();
    if $data[0] ~~ Numeric:D {
        for $data.List.kv -> $index, $value {
            die 'A response vector must contain only numeric values.' unless $value ~~ Numeric:D;
            @!data.push([$index, $value]);
        }
    }
    else {
        for $data.List -> $row {
            die 'Data rows must have at least two numeric values.'
            unless $row ~~ Positional:D && $row.elems >= 2
                    && $row[0] ~~ Numeric:D && $row[1] ~~ Numeric:D;
            @!data.push([$row[0], $row[1]]);
        }
    }

    return self
}

method take-data(--> Array:D) { @!data.Array }

method set-basis-functions($functions --> ::?CLASS:D) {
    die 'Basis functions must be a non-empty positional collection of callables.'
    unless $functions ~~ Positional:D && $functions.elems && $functions.all ~~ Callable:D;

    @!basis-functions = $functions.List;

    return self
}

method take-basis-functions(--> Array:D) { @!basis-functions.Array }

method set-regression-quantiles(%quantiles --> ::?CLASS:D) {
    die 'Regression quantiles must be a non-empty hash of callables.'
    unless %quantiles.elems && %quantiles.values.all ~~ Callable:D;

    %!regression-quantiles = %quantiles;

    return self
}

#| Take the regression quantile functions
method take-regression-quantiles(--> Hash:D) { %!regression-quantiles.Hash }

#| Set the pipeline value
method set-value($value --> ::?CLASS:D) { $!value = $value; self }

#| Take the pipeline value
method take-value() { $!value }

#| Echo the argument
method echo($message --> ::?CLASS:D) { say $message; self }

#==========================================================
# Summarization
#==========================================================

#| Data summary
method data-summary(--> ::?CLASS:D) {
    sub percentile(@values, Numeric:D $probability) {
        my @sorted = @values.sort;
        my $position = (@sorted.end * $probability).Num;
        my $lower = $position.floor;
        my $upper = $position.ceiling;
        @sorted[$lower] + ($position - $lower) * (@sorted[$upper] - @sorted[$lower])
    }
    my %summary;
    for ('regressor', 'value').kv -> $column, $name {
        my @values = @!data.map(*[$column]);
        %summary{$name} = {
            min => @values.min, '25%' => percentile(@values, .25),
            median => percentile(@values, .5), '75%' => percentile(@values, .75), max => @values.max,
        };
    }
    $!value = %summary;

    return self
}

#| Echo data summary
method echo-data-summary(Bool :$echo = True --> ::?CLASS:D) {
    self.data-summary;
    if $echo {
        say "Statistic    regressor | value";
        for <min 25% median 75% max> -> $statistic {
            say "$statistic\t" ~ %($!value)<regressor>{$statistic} ~ ' | ' ~ %($!value)<value>{$statistic};
        }
    }

    return self
}

#==========================================================
# Rescale
#==========================================================

#| Rescale the data
method rescale(Bool :$regressor = False, Bool :$value = True --> ::?CLASS:D) {
    for ($regressor, $value).kv -> $column, $enabled {
        next unless $enabled;
        my @values = @!data.map(*[$column]);
        my $minimum = @values.min;
        my $range = @values.max - $minimum;
        die 'Cannot rescale a constant data column.' if $range == 0;
        @!data[$_][$column] = (@!data[$_][$column] - $minimum) / $range for ^@!data.elems;
    }

    return self
}

#==========================================================
# Quantile regression
#==========================================================

#| Compute a Quantile Regression fit with a given basis of functions
method quantile-regression-fit($functions, :$probs, :$probabilities, *%options --> ::?CLASS:D) {
    self.set-basis-functions($functions);
    my @probs = self!probabilities($probabilities // $probs);
    my @fits = core-quantile-regression-fit(
        @!data, :basis(@!basis-functions), :probabilities(@probs), |%options
    );
    %!regression-quantiles = @probs Z=> @fits;
    $!value = %!regression-quantiles.Hash;

    return self
}

#| Compute a Quantile Regression fit with a B-spline basis
multi method quantile-regression($knots, $probabilities = Whatever, *%options --> ::?CLASS:D) {
    return self.quantile-regression(:$knots, :$probabilities, |%options)
}

#| Compute a Quantile Regression fit with a B-spline basis
multi method quantile-regression(:$knots, :probs(:$probabilities), :order(:$interpolation-order), *%options --> ::?CLASS:D) {
    my @probs = self!probabilities($probabilities);
    my $spline-order = $interpolation-order // 3;
    my @fits = core-quantile-regression(
        @!data, :$knots, :probabilities(@probs), :interpolation-order($spline-order), |%options
    );
    %!regression-quantiles = @probs Z=> @fits;
    $!value = %!regression-quantiles.Hash;
    return self
}

#==========================================================
# Least Squares
#==========================================================

method !solve(@matrix, @vector --> Array:D) {
    my @augmented = @matrix.kv.map(-> $row, @values { [|@values, @vector[$row]].Array }).Array;
    for ^@matrix.elems -> $pivot {
        my $best = $pivot;
        for ($pivot + 1 ..^ @augmented.elems) -> $candidate {
            $best = $candidate if abs(@augmented[$candidate][$pivot]) > abs(@augmented[$best][$pivot]);
        }
        die 'The least-squares basis is linearly dependent.' if abs(@augmented[$best][$pivot]) < 1e-12;
        (@augmented[$pivot], @augmented[$best]) = (@augmented[$best], @augmented[$pivot]);
        my $divisor = @augmented[$pivot][$pivot];
        @augmented[$pivot] »/=» $divisor;
        for ^@augmented.elems -> $row {
            next if $row == $pivot;
            my $factor = @augmented[$row][$pivot];
            @augmented[$row] »-=» (@augmented[$pivot] »*» $factor);
        }
    }
    return @augmented.map(*.tail).Array
}

#| Compute a Least Squares fit with a given basis of funcitons
method least-squares-fit($functions --> ::?CLASS:D) {
    self.set-basis-functions($functions);
    my $count = @!basis-functions.elems;
    my @normal = (0 xx $count).map({ (0 xx $count).Array }).Array;
    my @right = 0 xx $count;
    for @!data -> $row {
        my @basis = @!basis-functions.map(-> &function { function($row[0]) });
        die 'A basis function returned a non-numeric value.' unless @basis.all ~~ Numeric:D;
        for ^$count -> $i {
            @right[$i] += @basis[$i] * $row[1];
            for ^$count -> $j {
                @normal[$i][$j] += @basis[$i] * @basis[$j];
            }
        }
    }
    my @coefficients = self!solve(@normal, @right);
    %!regression-quantiles<mean> = -> $x { [+] @coefficients.kv.map(-> $i, $coefficient { $coefficient * @!basis-functions[$i]($x) }) };
    $!value = %!regression-quantiles.Hash;

    return self
}

#| Synonym of least-squares-fit
method fit($functions --> ::?CLASS:D) { self.least-squares-fit($functions) }

#==========================================================
# Workflow methods
#==========================================================

#| Evaluate the regression quantiles at given points
method evaluate($points --> ::?CLASS:D) {
    self!require-quantiles;
    my @points = $points ~~ Numeric:D ?? [$points,] !! $points.List;
    die 'Evaluation points must be numeric.' unless @points.all ~~ Numeric:D;
    my %values;
    for %!regression-quantiles.kv -> $probability, &fit {
        %values{$probability} = @points.map(&fit).Array;
    }
    $!value = %values;

    return self
}

#| Pick data points using regression quantiles as references
method pick-path-points(Numeric:D :$threshold, Bool :$pick-above-threshold = False, Bool :$relative-errors = False --> ::?CLASS:D) {
    die 'The threshold must be non-negative.' if $threshold < 0;
    self!require-quantiles;
    my %points;
    for %!regression-quantiles.kv -> $probability, &fit {
        %points{$probability} = @!data.grep(-> $row {
            my $error = abs($row[1] - fit($row[0]));
            $error /= abs($row[1]) if $relative-errors && $row[1] != 0;
            $pick-above-threshold ?? $error > $threshold !! $error <= $threshold
        }).Array;
    }
    $!value = %points;

    return self
}

#| Separate the data points using the regression quantiles
method separate(Bool :$cumulative = True, Bool :$fractions = False --> ::?CLASS:D) {
    self!require-quantiles;
    die 'Non-cumulative separation is not supported.' unless $cumulative;
    my %groups;
    for %!regression-quantiles.kv -> $probability, &fit {
        next if $probability eq 'mean';
        %groups{$probability} = @!data.grep(-> $row { $row[1] <= fit($row[0]) }).Array;
    }
    if $fractions {
        my %fractions;
        for %groups.kv -> $probability, @group {
            %fractions{$probability} = @group.elems / @!data.elems;
        }
        $!value = %fractions;
    }
    else {
        $!value = %groups;
    }

    return self
}

#| Compute fitting errors (residuals) per regression quantile
method errors(Bool :$relative-errors = False --> ::?CLASS:D) {
    self!require-quantiles;
    my %errors;
    for %!regression-quantiles.kv -> $probability, &fit {
        %errors{$probability} = @!data.map(-> $row {
            my $error = $row[1] - fit($row[0]);
            $error /= $row[1] if $relative-errors && $row[1] != 0;
            [$row[0], $error]
        }).Array;
    }
    $!value = %errors;

    return self
}

#| Find anomalies by residuals
method find-anomalies-by-residuals(Numeric :$threshold!, Bool :$relative-errors = False --> ::?CLASS:D) {
    self.pick-path-points(:$threshold, :pick-above-threshold, :$relative-errors);
    $!value = $!value.values.head.Array;

    return self
}

#| Find outliers
method outliers(--> ::?CLASS:D) {
    self!require-quantiles;
    my @probabilities = %!regression-quantiles.keys.grep(* ne 'mean').map(*.Num).sort;
    die 'At least two numeric regression quantiles are required.' unless @probabilities.elems >= 2;
    my (&bottom, &top) = %!regression-quantiles{@probabilities.head.Str}, %!regression-quantiles{@probabilities.tail.Str};
    $!value = {
        bottom => @!data.grep({ .[1] < bottom(.[0]) }).Array,
        top => @!data.grep({ .[1] > top(.[0]) }).Array,
    };

    return self
}

#| Compute conditional CDF
method conditional-cdf($points --> ::?CLASS:D) {
    self!require-quantiles;
    my @points = $points ~~ Numeric:D ?? [$points,] !! $points.List;

    my @probabilities = %!regression-quantiles.keys.grep(* ne 'mean').map(*.Num).sort;
    die 'At least two numeric regression quantiles are required.' unless @probabilities.elems >= 2;

    $!value = @points.map(-> $point {
        my @pairs = @probabilities.map({ [%!regression-quantiles{$_.Str}($point), $_] }).sort(*[0]);
        $point => sub ($value) {
            return @pairs.head[1] if $value <= @pairs.head[0];
            return @pairs.tail[1] if $value >= @pairs.tail[0];
            for ^@pairs.end -> $i {
                next unless @pairs[$i][0] <= $value <= @pairs[$i + 1][0];
                my ($x1, $p1) = @pairs[$i].List;
                my ($x2, $p2) = @pairs[$i + 1].List;
                return $x2 == $x1 ?? $p1 !! $p1 + ($value - $x1) * ($p2 - $p1) / ($x2 - $x1);
            }
        }
    }).Hash;

    return self
}

#| Simulate signal
proto method simulate($points) {*}

multi method simulate(Int:D $points) {
    die 'An integer point count must be larger than one.' unless $points > 1;

    # The stored data is always normalized to predictor-response pairs.
    my $minimum = @!data.map(*[0]).min;
    my $maximum = @!data.map(*[0]).max;
    my @simulation-points = (^$points).map({
        $minimum + ($maximum - $minimum) * $_ / ($points - 1)
    }).Array;
    self.simulate(@simulation-points)
}

multi method simulate(Positional:D $points --> ::?CLASS:D) {
    my @points = $points.List;
    die 'Simulation points must be a non-empty positional collection of numbers.'
    unless @points.elems && @points.all ~~ Numeric:D;

    # A least-squares fit can add the non-numeric key "mean".  
    # Simulation is defined by ordered regression quantiles only.
    my @probabilities = %!regression-quantiles.keys.grep(* ne 'mean').map(*.Num).sort;

    die 'Compute at least two numeric regression quantiles first.'
    unless @probabilities.elems >= 2;

    my @result;
    for @points -> $point {
        my @quantile-values = @probabilities.map(
            -> $probability { %!regression-quantiles{$probability.Str}($point) }
        ).Array;

        die 'Regression quantiles must evaluate to numeric values.'
        unless @quantile-values.all ~~ Numeric:D;

        # Select an adjacent quantile interval with probability proportional
        # to its width, then sample uniformly within that selected interval.
        my @weights = (^@quantile-values.end).map(
            -> $index { abs(@quantile-values[$index + 1] - @quantile-values[$index]) }
        ).Array;
        my $total-weight = [+] @weights;
        
        die 'Regression quantiles have zero total width at a simulation point.'
        unless $total-weight > 0;

        my $cutoff = rand * $total-weight;
        my $interval = @weights.end;
        my $cumulative-weight = 0e0;
        for @weights.kv -> $index, $weight {
            $cumulative-weight += $weight;
            if $cutoff < $cumulative-weight {
                $interval = $index;
                last;
            }
        }

        my $lower = @quantile-values[$interval];
        my $upper = @quantile-values[$interval + 1];
        @result.push([$point, $lower + rand * ($upper - $lower)]);
    }

    $!value = @result.Array;

    self
}


#==========================================================
# Plots
#==========================================================

#| Plot data and regression quantiles
method plot(
        :$title = 'Regression',      #= Plot title.
        Bool:D :$date-plot = False,  #= Whether to make a date plot or not.
        Numeric:D :$date-offset = 0, #= Use 2208988800 for epoch start at 1900-01-01 instead of 1970-01-01.
        :$data-color = 'Gray',       #= Color of the data points.
        Int:D :$width = 800,         #= Width of the plot.
        Int:D :$height = 600,        #= Height of the plot.
        :$format = 'html',           #= Format, one of 'html' or 'jupyter'.
        *%options,                   #= Additional options.
        --> ::?CLASS:D) {

    self!require-quantiles;

    my @names = 'x', 'data', |%!regression-quantiles.keys.sort;
    my @records = @!data.map(-> $row {
        my $x2 = $date-plot ?? DateTime.new(Instant.from-posix($row[0] - $date-offset)) !! $row[0];
        (@names Z=> [$x2, $row[1], |%!regression-quantiles.keys.sort.map({ %!regression-quantiles{$_}($row[0]) })]).Hash
    }).Array;

    $!value = js-google-charts('Combo',
            @records,
            column-names => @names,
            :$title,
            :$width, :$height,
            :$format,
            seriesType => 'line',
            series => { 0 => { type => 'scatter', color => $data-color, pointSize => 3, lineWidth => 0 } },
            |%options);

    return self
}

#| Plot data and regression quantiles using date horizontal axis
method date-list-plot(*%args) {
    return self.plot(:date-plot, |%args)
}

#| Outliers plot
method outliers-plot(
        :$title = 'Regression',      #= Plot title.
        Bool:D :$date-plot = False,  #= Whether to make a date plot or not.
        Numeric:D :$date-offset = 0, #= Use 2208988800 for epoch start at 1900-01-01 instead of 1970-01-01.
        :$data-color = 'Gray',       #= Color of the data points.
        Int:D :$width = 800,         #= Width of the plot.
        Int:D :$height = 600,        #= Height of the plot.
        :$format = 'html',           #= Format, one of 'html' or 'jupyter'.
        *%options,                   #= Additional options.
                     ) {

    my &bottom-func = %!regression-quantiles.grep(*.key ne 'mean').sort(*.key.Num).head.value;
    my &top-func = %!regression-quantiles.grep(*.key ne 'mean').sort(*.key.Num).tail.value;
    # Get the bottom and top outliers
    my @wide-data = @!data.kv.map(-> $i, @point {
        %(
            x => $date-plot ?? DateTime.new(Instant.from-posix(@point.head - $date-offset)) !! @point.head,
            data => @point.tail,
            bottom => bottom-func(@point.head),
            top => top-func(@point.head),
            bottom-outlier => bottom-func(@point.head) > @point.tail ?? @point.tail !! Nil,
            top-outlier => top-func(@point.head) < @point.tail ?? @point.tail !! Nil,
        )
    }).Array;

    my @field-names = <x data bottom top bottom-outlier top-outlier>;
    my @column-names = ['x', 'data', 'bottom regression quantile', 'top regression quantile', 'bottom outlier', 'top outlier'];

    # Plot
    $!value = js-google-charts('Scatter',
            @wide-data.map(*{|@field-names}),
            :@column-names,
            :$title,
            series => {
                0 => {color => $data-color, pointSize => 2, title => 'data'},
                1 => {color => 'LightBlue', type => 'line', lineWidth => 2, pointSize => 0},
                2 => {color => 'Pink', type => 'line', lineWidth => 2, pointSize => 0},
                3 => {color => 'DodgerBlue', pointSize => 3},
                4 => {color => 'Red', pointSize => 3},
            },
            legend => {
                position => 'right',
                textStyle => {
                    color => 'Gray',
                    fontSize => 12
                }
            },
            :$format,
            :$width,
            :$height,
            |%options);

    return self
}

#| Errors plot
method error-plots(
        :$title = 'Regression',              #= Plot title.
        Str:D :type($plot-type) = 'Scatter', #= Plot type; one of "Scatter" or 'Line'.
        Bool:D :$date-plot = False,          #= Whether to make a date plot or not.
        Numeric:D :$date-offset = 0,         #= Use 2208988800 for epoch start at 1900-01-01 instead of 1970-01-01.
        Int:D :$width = 800,                 #= Width of the plot.
        Int:D :$height = 600,                #= Height of the plot.
        :$format = 'html',                   #= Format, one of 'html' or 'jupyter'.
        *%options,                           #= Additional options.
                    ) {
    # Check plot type
    die 'The value $plot-type is expected to be one "Line" or "Scatter".' unless $plot-type ∈ <Line Scatter>;

    # Compute errors
    self.errors;

    # Modify options
    my %opts = %options.grep(*.key ∉ <title div-id>);

    # In order to handle plot, a copy of the errors is made.
    my %errors = self.take-value.clone;

    # Modify x-axis
    if $date-plot {
        %errors = %errors.kv.map( -> $prob, @points { $prob => @points.map({ [DateTime.new(Instant.from-posix($_.head - $date-offset)), $_.tail] }).Array });
    }

    # Plots
    my @res =
            %errors.sort(*.key).map({
                js-google-charts($plot-type,
                        $_.value,
                        title => (%options<title> // "Errors of the regression quantile at") ~ " " ~ $_.key,
                        column-names => ['x', $_.key],
                        :$width,
                        :$height,
                        :$format,
                        |%opts,
                        div-id => "{%options<div-id> // 'errors'}-{ $_.key }")
            });

    $!value = @res;

    return self
}

#| Conditional CDF plot(s)
multi method conditional-cdf-plots($point, *%args) {
    return self.conditional-cdf-plot([$point, ], |%args)
}

#| Conditional CDF plot(s)
multi method conditional-cdf-plots(
        @points,
        :$title = '',              #= Plot title.
        Int:D :$width = 800,       #= Width of the plot.
        Int:D :$height = 600,      #= Height of the plot.
        :$format = 'html',         #= Format, one of 'html' or 'jupyter'.
        *%options,                 #= Additional options.
                    ) {
    die 'A number or a list of numbers is expected as first argument.' unless @points.all ~~ Numeric:D;

    my @probabilities = %!regression-quantiles.keys.grep(* ne 'mean').map(*.Num).sort;

    die 'Compute at least two numeric regression quantiles first.'
    unless @probabilities.elems >= 2;

    my @q-funcs = %!regression-quantiles{|@probabilities};

    # For each point compute the its conditional CDF
    my %cdfs = do for @points -> $p { 
        my @q-func-values = @q-funcs.kv.map(-> $i, &f { [&f($p), @probabilities[$i]] });
        $p => @q-func-values
    }

    # Modify options
    my %opts = %options.grep(*.key ∉ <title div-id>);

    # Plots
    my @res =
            %cdfs.map({
                js-google-charts('Line',
                        $_.value,
                        title => (%options<title> // "Conditional CDF at") ~ " " ~ $_.key,
                        column-names => ['x', $_.key],
                        :$width,
                        :$height,
                        :$format,
                        |%opts,
                        div-id => "{%options<div-id> // 'cdfs'}-{ $_.key }")
            });

    $!value = @res;

    return self
}

#==========================================================
# Representation
#==========================================================

method Str(::?CLASS:D:-->Str) {
    self.gist
}

multi method Hash(::?CLASS:D:-->Map:D) {
    my %h =
            :@!data,
            :@!basis-functions,
            :%!regression-quantiles,
            :$!value
            ;
    return %h
}

multi method gist(::?CLASS:D:-->Str) {
    'Math::Regressionizer(' ~
            ' data => ' ~ @!data.elems ~ ' points'~
            ', basis-functions => ' ~ @!basis-functions.elems ~ ' count'~
            ', regression-quantiles => ' ~ (%!regression-quantiles ?? 'at (' ~ (%!regression-quantiles.keys.sort).join(', ') ~ ')' !! '()') ~
            ', value => ' ~ ($!value.defined ?? $!value.WHAT.raku !! $!value.raku) ~
            ')'
}
