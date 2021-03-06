# This file was automatically generated by Dist::Zilla::Plugin::MakeMaker v6.009.
use strict;
use warnings;

use 5.010000;

use ExtUtils::MakeMaker;

my %WriteMakefileArgs = (
  "ABSTRACT" => "Run Tassel's version 2 of the GBS pipeline using parameters from a json configuration file",
  "AUTHOR" => "Christopher Bottoms <molecules <at> cpan <dot> org>",
  "BUILD_REQUIRES" => {
    "Module::Build" => "0.28"
  },
  "CONFIGURE_REQUIRES" => {
    "ExtUtils::MakeMaker" => 0,
    "Module::Build" => "0.28"
  },
  "DISTNAME" => "GBSV2-Tassel-Pipeline",
  "LICENSE" => "perl",
  "MIN_PERL_VERSION" => "5.010000",
  "NAME" => "GBSV2::Tassel::Pipeline",
  "PREREQ_PM" => {
    "Cwd" => 0,
    "Exporter" => 0,
    "File::Glob" => 0,
    "File::Slurper" => 0,
    "JSON" => 0,
    "autodie" => 0,
    "strict" => 0,
    "warnings" => 0
  },
  "TEST_REQUIRES" => {
    "Test2::Bundle::Extended" => 0,
    "Test::LongString" => 0,
    "lib" => 0
  },
  "VERSION" => "0.0009",
  "test" => {
    "TESTS" => "t/*.t"
  }
);


my %FallbackPrereqs = (
  "Cwd" => 0,
  "Exporter" => 0,
  "File::Glob" => 0,
  "File::Slurper" => 0,
  "JSON" => 0,
  "Module::Build" => "0.28",
  "Test2::Bundle::Extended" => 0,
  "Test::LongString" => 0,
  "autodie" => 0,
  "lib" => 0,
  "strict" => 0,
  "warnings" => 0
);


unless ( eval { ExtUtils::MakeMaker->VERSION(6.63_03) } ) {
  delete $WriteMakefileArgs{TEST_REQUIRES};
  delete $WriteMakefileArgs{BUILD_REQUIRES};
  $WriteMakefileArgs{PREREQ_PM} = \%FallbackPrereqs;
}

delete $WriteMakefileArgs{CONFIGURE_REQUIRES}
  unless eval { ExtUtils::MakeMaker->VERSION(6.52) };

WriteMakefile(%WriteMakefileArgs);
