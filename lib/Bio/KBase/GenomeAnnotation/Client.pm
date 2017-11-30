package Bio::KBase::GenomeAnnotation::Client;

use POSIX;
use strict;
use Data::Dumper;
use URI;

my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use P3AuthToken;



sub new
{
    my($class, $url, @args) = @_;
    
    if (!defined($url))
    {
	$url = 'https://p3.theseed.org/services/genome_annotation';
    }

    my $self = {
	client => Bio::KBase::GenomeAnnotation::Client::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = P3AuthToken->new(@args);
	
	if (my $token_str = $token->token())
	{
	    $self->{token} = $token_str;
	    $self->{client}->{token} = $token_str;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    $ua->agent("Bio::KBase::GenomeAnnotation::Client UserAgent");
    bless $self, $class;
    return $self;
}




sub genome_ids_to_genomes
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function genome_ids_to_genomes (received $n, expecting 1)";
    }
    {
	my($ids) = @args;

	my @_bad_arguments;
        (ref($ids) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 1 \"ids\" (value was \"$ids\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genome_ids_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.genome_ids_to_genomes",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking genome_ids_to_genomes:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method genome_ids_to_genomes: " .  $self->{client}->status_line;
    }
}



sub create_genome
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function create_genome (received $n, expecting 1)";
    }
    {
	my($metadata) = @args;

	my @_bad_arguments;
        (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"metadata\" (value was \"$metadata\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to create_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.create_genome",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking create_genome:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method create_genome: " .  $self->{client}->status_line;
    }
}



sub create_genome_from_genbank
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function create_genome_from_genbank (received $n, expecting 1)";
    }
    {
	my($gb_data) = @args;

	my @_bad_arguments;
        (!ref($gb_data)) or push(@_bad_arguments, "Invalid type for argument 1 \"gb_data\" (value was \"$gb_data\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to create_genome_from_genbank:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.create_genome_from_genbank",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking create_genome_from_genbank:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method create_genome_from_genbank: " .  $self->{client}->status_line;
    }
}



sub create_genome_from_SEED
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function create_genome_from_SEED (received $n, expecting 1)";
    }
    {
	my($genome_id) = @args;

	my @_bad_arguments;
        (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument 1 \"genome_id\" (value was \"$genome_id\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to create_genome_from_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.create_genome_from_SEED",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking create_genome_from_SEED:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method create_genome_from_SEED: " .  $self->{client}->status_line;
    }
}



sub create_genome_from_RAST
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function create_genome_from_RAST (received $n, expecting 1)";
    }
    {
	my($genome_or_job_id) = @args;

	my @_bad_arguments;
        (!ref($genome_or_job_id)) or push(@_bad_arguments, "Invalid type for argument 1 \"genome_or_job_id\" (value was \"$genome_or_job_id\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to create_genome_from_RAST:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.create_genome_from_RAST",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking create_genome_from_RAST:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method create_genome_from_RAST: " .  $self->{client}->status_line;
    }
}



sub set_metadata
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function set_metadata (received $n, expecting 2)";
    }
    {
	my($genome_in, $metadata) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"metadata\" (value was \"$metadata\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.set_metadata",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking set_metadata:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method set_metadata: " .  $self->{client}->status_line;
    }
}



sub add_contigs
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function add_contigs (received $n, expecting 2)";
    }
    {
	my($genome_in, $contigs) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"contigs\" (value was \"$contigs\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.add_contigs",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking add_contigs:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method add_contigs: " .  $self->{client}->status_line;
    }
}



sub add_contigs_from_handle
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function add_contigs_from_handle (received $n, expecting 2)";
    }
    {
	my($genome_in, $contigs) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"contigs\" (value was \"$contigs\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to add_contigs_from_handle:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.add_contigs_from_handle",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking add_contigs_from_handle:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method add_contigs_from_handle: " .  $self->{client}->status_line;
    }
}



sub add_features
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function add_features (received $n, expecting 2)";
    }
    {
	my($genome_in, $features) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($features) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"features\" (value was \"$features\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to add_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.add_features",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking add_features:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method add_features: " .  $self->{client}->status_line;
    }
}



sub genomeTO_to_reconstructionTO
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function genomeTO_to_reconstructionTO (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.genomeTO_to_reconstructionTO",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking genomeTO_to_reconstructionTO:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method genomeTO_to_reconstructionTO: " .  $self->{client}->status_line;
    }
}



sub genomeTO_to_feature_data
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function genomeTO_to_feature_data (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.genomeTO_to_feature_data",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking genomeTO_to_feature_data:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method genomeTO_to_feature_data: " .  $self->{client}->status_line;
    }
}



sub reconstructionTO_to_roles
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function reconstructionTO_to_roles (received $n, expecting 1)";
    }
    {
	my($reconstructionTO) = @args;

	my @_bad_arguments;
        (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"reconstructionTO\" (value was \"$reconstructionTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.reconstructionTO_to_roles",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking reconstructionTO_to_roles:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method reconstructionTO_to_roles: " .  $self->{client}->status_line;
    }
}



sub reconstructionTO_to_subsystems
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function reconstructionTO_to_subsystems (received $n, expecting 1)";
    }
    {
	my($reconstructionTO) = @args;

	my @_bad_arguments;
        (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"reconstructionTO\" (value was \"$reconstructionTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.reconstructionTO_to_subsystems",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking reconstructionTO_to_subsystems:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method reconstructionTO_to_subsystems: " .  $self->{client}->status_line;
    }
}



sub assign_functions_to_CDSs
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function assign_functions_to_CDSs (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.assign_functions_to_CDSs",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking assign_functions_to_CDSs:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method assign_functions_to_CDSs: " .  $self->{client}->status_line;
    }
}



sub annotate_genome
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_genome (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_genome",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_genome:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_genome: " .  $self->{client}->status_line;
    }
}



sub call_selenoproteins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_selenoproteins (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_selenoproteins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_selenoproteins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_selenoproteins: " .  $self->{client}->status_line;
    }
}



sub call_pyrrolysoproteins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_pyrrolysoproteins (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_pyrrolysoproteins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_pyrrolysoproteins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_pyrrolysoproteins: " .  $self->{client}->status_line;
    }
}



sub call_features_selenoprotein
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_selenoprotein (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_selenoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_selenoprotein",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_selenoprotein:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_selenoprotein: " .  $self->{client}->status_line;
    }
}



sub call_features_pyrrolysoprotein
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_pyrrolysoprotein (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_pyrrolysoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_pyrrolysoprotein",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_pyrrolysoprotein:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_pyrrolysoprotein: " .  $self->{client}->status_line;
    }
}



sub call_features_insertion_sequences
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_insertion_sequences (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_insertion_sequences:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_insertion_sequences",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_insertion_sequences:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_insertion_sequences: " .  $self->{client}->status_line;
    }
}



sub call_features_rRNA_SEED
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_rRNA_SEED (received $n, expecting 2)";
    }
    {
	my($genome_in, $types) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"types\" (value was \"$types\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_rRNA_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_rRNA_SEED",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_rRNA_SEED:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_rRNA_SEED: " .  $self->{client}->status_line;
    }
}



sub call_features_tRNA_trnascan
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_tRNA_trnascan (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_tRNA_trnascan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_tRNA_trnascan",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_tRNA_trnascan:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_tRNA_trnascan: " .  $self->{client}->status_line;
    }
}



sub call_RNAs
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_RNAs (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_RNAs",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_RNAs:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_RNAs: " .  $self->{client}->status_line;
    }
}



sub call_features_CDS_glimmer3
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_CDS_glimmer3 (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_CDS_glimmer3",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_CDS_glimmer3:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_CDS_glimmer3: " .  $self->{client}->status_line;
    }
}



sub call_features_CDS_prodigal
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_CDS_prodigal (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_CDS_prodigal:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_CDS_prodigal",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_CDS_prodigal:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_CDS_prodigal: " .  $self->{client}->status_line;
    }
}



sub call_features_CDS_genemark
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_CDS_genemark (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_CDS_genemark:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_CDS_genemark",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_CDS_genemark:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_CDS_genemark: " .  $self->{client}->status_line;
    }
}



sub call_features_CDS_SEED_projection
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_CDS_SEED_projection (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_CDS_SEED_projection:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_CDS_SEED_projection",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_CDS_SEED_projection:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_CDS_SEED_projection: " .  $self->{client}->status_line;
    }
}



sub call_features_CDS_FragGeneScan
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_CDS_FragGeneScan (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_CDS_FragGeneScan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_CDS_FragGeneScan",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_CDS_FragGeneScan:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_CDS_FragGeneScan: " .  $self->{client}->status_line;
    }
}



sub call_features_repeat_region_SEED
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_repeat_region_SEED (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_repeat_region_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_repeat_region_SEED",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_repeat_region_SEED:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_repeat_region_SEED: " .  $self->{client}->status_line;
    }
}



sub call_features_prophage_phispy
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_prophage_phispy (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_prophage_phispy:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_prophage_phispy",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_prophage_phispy:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_prophage_phispy: " .  $self->{client}->status_line;
    }
}



sub call_features_scan_for_matches
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 3)
    {
        die "Invalid argument count for function call_features_scan_for_matches (received $n, expecting 3)";
    }
    {
	my($genome_in, $pattern, $feature_type) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (!ref($pattern)) or push(@_bad_arguments, "Invalid type for argument 2 \"pattern\" (value was \"$pattern\")");
        (!ref($feature_type)) or push(@_bad_arguments, "Invalid type for argument 3 \"feature_type\" (value was \"$feature_type\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_scan_for_matches:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_scan_for_matches",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_scan_for_matches:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_scan_for_matches: " .  $self->{client}->status_line;
    }
}



sub call_features_assembly_gap
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_assembly_gap (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_assembly_gap:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_assembly_gap",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_assembly_gap:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_assembly_gap: " .  $self->{client}->status_line;
    }
}



sub annotate_proteins_similarity
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function annotate_proteins_similarity (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins_similarity:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_proteins_similarity",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_proteins_similarity:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_proteins_similarity: " .  $self->{client}->status_line;
    }
}



sub annotate_proteins_phage
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function annotate_proteins_phage (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins_phage:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_proteins_phage",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_proteins_phage:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_proteins_phage: " .  $self->{client}->status_line;
    }
}



sub annotate_proteins_kmer_v1
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function annotate_proteins_kmer_v1 (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_proteins_kmer_v1",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_proteins_kmer_v1:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_proteins_kmer_v1: " .  $self->{client}->status_line;
    }
}



sub annotate_proteins_kmer_v2
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function annotate_proteins_kmer_v2 (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_proteins_kmer_v2",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_proteins_kmer_v2:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_proteins_kmer_v2: " .  $self->{client}->status_line;
    }
}



sub resolve_overlapping_features
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function resolve_overlapping_features (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to resolve_overlapping_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.resolve_overlapping_features",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking resolve_overlapping_features:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method resolve_overlapping_features: " .  $self->{client}->status_line;
    }
}



sub propagate_genbank_feature_metadata
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function propagate_genbank_feature_metadata (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to propagate_genbank_feature_metadata:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.propagate_genbank_feature_metadata",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking propagate_genbank_feature_metadata:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method propagate_genbank_feature_metadata: " .  $self->{client}->status_line;
    }
}



sub call_features_ProtoCDS_kmer_v1
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_ProtoCDS_kmer_v1 (received $n, expecting 2)";
    }
    {
	my($genomeTO, $params) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_ProtoCDS_kmer_v1",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_ProtoCDS_kmer_v1:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_ProtoCDS_kmer_v1: " .  $self->{client}->status_line;
    }
}



sub call_features_ProtoCDS_kmer_v2
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function call_features_ProtoCDS_kmer_v2 (received $n, expecting 2)";
    }
    {
	my($genome_in, $params) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_ProtoCDS_kmer_v2",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_ProtoCDS_kmer_v2:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_ProtoCDS_kmer_v2: " .  $self->{client}->status_line;
    }
}



sub enumerate_special_protein_databases
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function enumerate_special_protein_databases (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.enumerate_special_protein_databases",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking enumerate_special_protein_databases:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method enumerate_special_protein_databases: " .  $self->{client}->status_line;
    }
}



sub compute_special_proteins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function compute_special_proteins (received $n, expecting 2)";
    }
    {
	my($genome_in, $database_names) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($database_names) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"database_names\" (value was \"$database_names\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to compute_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.compute_special_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking compute_special_proteins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method compute_special_proteins: " .  $self->{client}->status_line;
    }
}



sub annotate_special_proteins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_special_proteins (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_special_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_special_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_special_proteins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_special_proteins: " .  $self->{client}->status_line;
    }
}



sub annotate_families_figfam_v1
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_families_figfam_v1 (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_families_figfam_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_families_figfam_v1",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_families_figfam_v1:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_families_figfam_v1: " .  $self->{client}->status_line;
    }
}



sub annotate_families_patric
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_families_patric (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_families_patric:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_families_patric",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_families_patric:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_families_patric: " .  $self->{client}->status_line;
    }
}



sub annotate_null_to_hypothetical
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_null_to_hypothetical (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_null_to_hypothetical:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_null_to_hypothetical",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_null_to_hypothetical:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_null_to_hypothetical: " .  $self->{client}->status_line;
    }
}



sub annotate_strain_type_MLST
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_strain_type_MLST (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_strain_type_MLST:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_strain_type_MLST",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_strain_type_MLST:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_strain_type_MLST: " .  $self->{client}->status_line;
    }
}



sub compute_cdd
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function compute_cdd (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to compute_cdd:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.compute_cdd",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking compute_cdd:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method compute_cdd: " .  $self->{client}->status_line;
    }
}



sub annotate_proteins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function annotate_proteins (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.annotate_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking annotate_proteins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method annotate_proteins: " .  $self->{client}->status_line;
    }
}



sub estimate_crude_phylogenetic_position_kmer
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function estimate_crude_phylogenetic_position_kmer (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to estimate_crude_phylogenetic_position_kmer:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.estimate_crude_phylogenetic_position_kmer",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking estimate_crude_phylogenetic_position_kmer:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method estimate_crude_phylogenetic_position_kmer: " .  $self->{client}->status_line;
    }
}



sub find_close_neighbors
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function find_close_neighbors (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.find_close_neighbors",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking find_close_neighbors:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method find_close_neighbors: " .  $self->{client}->status_line;
    }
}



sub call_features_strep_suis_repeat
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_strep_suis_repeat (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_strep_suis_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_strep_suis_repeat",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_strep_suis_repeat:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_strep_suis_repeat: " .  $self->{client}->status_line;
    }
}



sub call_features_strep_pneumo_repeat
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_strep_pneumo_repeat (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_strep_pneumo_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_strep_pneumo_repeat",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_strep_pneumo_repeat:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_strep_pneumo_repeat: " .  $self->{client}->status_line;
    }
}



sub call_features_crispr
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function call_features_crispr (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.call_features_crispr",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking call_features_crispr:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method call_features_crispr: " .  $self->{client}->status_line;
    }
}



sub update_functions
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 3)
    {
        die "Invalid argument count for function update_functions (received $n, expecting 3)";
    }
    {
	my($genome_in, $functions, $event) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($functions) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"functions\" (value was \"$functions\")");
        (ref($event) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 3 \"event\" (value was \"$event\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to update_functions:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.update_functions",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking update_functions:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method update_functions: " .  $self->{client}->status_line;
    }
}



sub renumber_features
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function renumber_features (received $n, expecting 1)";
    }
    {
	my($genome_in) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to renumber_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.renumber_features",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking renumber_features:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method renumber_features: " .  $self->{client}->status_line;
    }
}



sub classify_amr
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function classify_amr (received $n, expecting 1)";
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to classify_amr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.classify_amr",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking classify_amr:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method classify_amr: " .  $self->{client}->status_line;
    }
}



sub export_genome
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 3)
    {
        die "Invalid argument count for function export_genome (received $n, expecting 3)";
    }
    {
	my($genome_in, $format, $feature_types) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (!ref($format)) or push(@_bad_arguments, "Invalid type for argument 2 \"format\" (value was \"$format\")");
        (ref($feature_types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 3 \"feature_types\" (value was \"$feature_types\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to export_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.export_genome",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking export_genome:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method export_genome: " .  $self->{client}->status_line;
    }
}



sub enumerate_classifiers
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function enumerate_classifiers (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.enumerate_classifiers",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking enumerate_classifiers:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method enumerate_classifiers: " .  $self->{client}->status_line;
    }
}



sub query_classifier_groups
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function query_classifier_groups (received $n, expecting 1)";
    }
    {
	my($classifier) = @args;

	my @_bad_arguments;
        (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument 1 \"classifier\" (value was \"$classifier\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to query_classifier_groups:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.query_classifier_groups",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking query_classifier_groups:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method query_classifier_groups: " .  $self->{client}->status_line;
    }
}



sub query_classifier_taxonomies
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function query_classifier_taxonomies (received $n, expecting 1)";
    }
    {
	my($classifier) = @args;

	my @_bad_arguments;
        (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument 1 \"classifier\" (value was \"$classifier\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to query_classifier_taxonomies:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.query_classifier_taxonomies",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking query_classifier_taxonomies:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method query_classifier_taxonomies: " .  $self->{client}->status_line;
    }
}



sub classify_into_bins
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function classify_into_bins (received $n, expecting 2)";
    }
    {
	my($classifier, $dna_input) = @args;

	my @_bad_arguments;
        (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument 1 \"classifier\" (value was \"$classifier\")");
        (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"dna_input\" (value was \"$dna_input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to classify_into_bins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.classify_into_bins",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking classify_into_bins:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method classify_into_bins: " .  $self->{client}->status_line;
    }
}



sub classify_full
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function classify_full (received $n, expecting 2)";
    }
    {
	my($classifier, $dna_input) = @args;

	my @_bad_arguments;
        (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument 1 \"classifier\" (value was \"$classifier\")");
        (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 2 \"dna_input\" (value was \"$dna_input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to classify_full:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.classify_full",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking classify_full:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method classify_full: " .  $self->{client}->status_line;
    }
}



sub default_workflow
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function default_workflow (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.default_workflow",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking default_workflow:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method default_workflow: " .  $self->{client}->status_line;
    }
}



sub enumerate_workflows
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function enumerate_workflows (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.enumerate_workflows",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking enumerate_workflows:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method enumerate_workflows: " .  $self->{client}->status_line;
    }
}



sub complete_workflow_template
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function complete_workflow_template (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.complete_workflow_template",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking complete_workflow_template:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method complete_workflow_template: " .  $self->{client}->status_line;
    }
}



sub run_pipeline
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function run_pipeline (received $n, expecting 2)";
    }
    {
	my($genome_in, $workflow) = @args;

	my @_bad_arguments;
        (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genome_in\" (value was \"$genome_in\")");
        (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"workflow\" (value was \"$workflow\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.run_pipeline",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking run_pipeline:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method run_pipeline: " .  $self->{client}->status_line;
    }
}



sub pipeline_batch_start
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 2)
    {
        die "Invalid argument count for function pipeline_batch_start (received $n, expecting 2)";
    }
    {
	my($genomes, $workflow) = @args;

	my @_bad_arguments;
        (ref($genomes) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument 1 \"genomes\" (value was \"$genomes\")");
        (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 2 \"workflow\" (value was \"$workflow\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to pipeline_batch_start:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.pipeline_batch_start",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking pipeline_batch_start:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method pipeline_batch_start: " .  $self->{client}->status_line;
    }
}



sub pipeline_batch_status
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
        die "Invalid argument count for function pipeline_batch_status (received $n, expecting 1)";
    }
    {
	my($batch_id) = @args;

	my @_bad_arguments;
        (!ref($batch_id)) or push(@_bad_arguments, "Invalid type for argument 1 \"batch_id\" (value was \"$batch_id\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to pipeline_batch_status:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    die $msg;
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.pipeline_batch_status",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking pipeline_batch_status:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method pipeline_batch_status: " .  $self->{client}->status_line;
    }
}



sub pipeline_batch_enumerate_batches
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 0)
    {
        die "Invalid argument count for function pipeline_batch_enumerate_batches (received $n, expecting 0)";
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "GenomeAnnotation.pipeline_batch_enumerate_batches",
	params => \@args,
    });
    if ($result) {
	if ($result->{error}) {
	    die "Error $result->{error}->{code} invoking pipeline_batch_enumerate_batches:\n" . $self->{client}->json->encode($result->{error}->{error});
	} else {
	    return wantarray ? @{$result->{result}} : $result->{result}->[0];
	}
    } else {
	die "Error invoking method pipeline_batch_enumerate_batches: " .  $self->{client}->status_line;
    }
}




package Bio::KBase::GenomeAnnotation::Client::RpcClient;
use POSIX;
use strict;
use LWP::UserAgent;
use JSON::XS;

BEGIN {
    for my $method (qw/uri ua json content_type version id allow_call status_line/) {
	eval qq|
	    sub $method {
		\$_[0]->{$method} = \$_[1] if defined \$_[1];
		\$_[0]->{$method};
	    }
	    |;
	}
    }

sub new
{
    my($class) = @_;

    my $ua = LWP::UserAgent->new();
    my $json = JSON::XS->new->allow_nonref->utf8;
    
    my $self = {
	ua => $ua,
	json => $json,
    };
    return bless $self, $class;
}

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    my @retries = (1, 2, 5, 10, 20, 60, 60, 60, 60, 60, 60);
    my %codes_to_retry =  map { $_ => 1 } qw(110 408 502 503 504 200) ;
    my $n_retries;

    while (1)

    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

	#
	# Bail early on success.
	#
	if ($result->is_success)
	{
	    if ($n_retries)
	    {
		print STDERR strftime("%F %T", localtime), ": Request succeeded after $n_retries retries\n";
	    }
	    last;
	}
	$n_retries++;

	#
	# Failure. See if we need to retry and loop, or bail with
	# a permanent failure.
	#
	
        my $code = $result->code;
	my $msg = $result->message;
	my $want_retry = 0;
	if ($codes_to_retry{$code})
	{
	    $want_retry = 1;
	}
	elsif ($code eq 500 && defined( $result->header('client-warning') )
	       && $result->header('client-warning') eq 'Internal response')
	{
	    #
	    # Handle errors that were not thrown by the web
	    # server but rather picked up by the client library.
	    #
	    # If we got a client timeout or connection refused, let us retry.
	    #
	    
	    if ($msg =~ /timeout|connection refused/i)
	    {
		$want_retry = 1;
	    }
	    
	}
	
        if (!$want_retry || @retries == 0) {
	    last;
        }
	
        #
        # otherwise, sleep & loop.
        #
        my $retry_time = shift(@retries);
        print STDERR strftime("%F %T", localtime), ": Request failed with code=$code msg=$msg, sleeping $retry_time and retrying\n";
        sleep($retry_time);

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success || $result->content_type eq 'application/json') {

	my $txt = $result->content;

        return unless($txt); # notification?

	my $obj = eval { $self->json->decode($txt); };

	if (!$obj)
	{
	    die "Error parsing result: $@";
	}

	return $obj;
    }
    else {
        return;
    }
}

sub _get {
    my ($self, $uri) = @_;
    $self->ua->get(
		   $uri,
		   Accept         => 'application/json',
		  );
}

sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Legacy::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
