// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.





#include "variant_filter_parser.h"


// Unqualified %code blocks.
#line 36 "variant_filter_parser.ypp"

#include "variant_filter_driver.h"

#line 50 "variant_filter_parser.cpp"


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif


// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yy_stack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 15 "variant_filter_parser.ypp"
namespace kim {
#line 143 "variant_filter_parser.cpp"

  /// Build a parser object.
  VariantFilterParser::VariantFilterParser (VariantFilterDriver &driver_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      yy_lac_established_ (false),
      driver (driver_yyarg)
  {}

  VariantFilterParser::~VariantFilterParser ()
  {}

  VariantFilterParser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------.
  | symbol.  |
  `---------*/



  // by_state.
  VariantFilterParser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  VariantFilterParser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  VariantFilterParser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  VariantFilterParser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  VariantFilterParser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  VariantFilterParser::symbol_kind_type
  VariantFilterParser::by_state::kind () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST (symbol_kind_type, yystos_[+state]);
  }

  VariantFilterParser::stack_symbol_type::stack_symbol_type ()
  {}

  VariantFilterParser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.YY_MOVE_OR_COPY< bool > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.YY_MOVE_OR_COPY< double > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.YY_MOVE_OR_COPY< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.YY_MOVE_OR_COPY< kim::VariantFilter::BooleanOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.YY_MOVE_OR_COPY< kim::VariantFilter::ComparisonOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.YY_MOVE_OR_COPY< kim::VariantFilter::ExpressionOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.YY_MOVE_OR_COPY< kim::VariantFilter::FilterOperation > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_num_op: // num_op
        value.YY_MOVE_OR_COPY< kim::VariantFilter::NumericalOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_str_op: // str_op
        value.YY_MOVE_OR_COPY< kim::VariantFilter::StringOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.YY_MOVE_OR_COPY< std::string > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  VariantFilterParser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.move< bool > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.move< double > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.move< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.move< kim::VariantFilter::BooleanOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.move< kim::VariantFilter::ComparisonOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.move< kim::VariantFilter::ExpressionOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.move< kim::VariantFilter::FilterOperation > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_num_op: // num_op
        value.move< kim::VariantFilter::NumericalOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_str_op: // str_op
        value.move< kim::VariantFilter::StringOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.move< std::string > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  VariantFilterParser::stack_symbol_type&
  VariantFilterParser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.copy< bool > (that.value);
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.copy< double > (that.value);
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.copy< int > (that.value);
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.copy< kim::VariantFilter::BooleanOperator > (that.value);
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.copy< kim::VariantFilter::ComparisonOperator > (that.value);
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.copy< kim::VariantFilter::ExpressionOperator > (that.value);
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.copy< kim::VariantFilter::FilterOperation > (that.value);
        break;

      case symbol_kind::S_num_op: // num_op
        value.copy< kim::VariantFilter::NumericalOperator > (that.value);
        break;

      case symbol_kind::S_str_op: // str_op
        value.copy< kim::VariantFilter::StringOperator > (that.value);
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.copy< std::string > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }

  VariantFilterParser::stack_symbol_type&
  VariantFilterParser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.move< bool > (that.value);
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.move< double > (that.value);
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.move< int > (that.value);
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.move< kim::VariantFilter::BooleanOperator > (that.value);
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.move< kim::VariantFilter::ComparisonOperator > (that.value);
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.move< kim::VariantFilter::ExpressionOperator > (that.value);
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.move< kim::VariantFilter::FilterOperation > (that.value);
        break;

      case symbol_kind::S_num_op: // num_op
        value.move< kim::VariantFilter::NumericalOperator > (that.value);
        break;

      case symbol_kind::S_str_op: // str_op
        value.move< kim::VariantFilter::StringOperator > (that.value);
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.move< std::string > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  VariantFilterParser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  VariantFilterParser::yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YY_USE (yyoutput);
    if (yysym.empty ())
      yyo << "empty symbol";
    else
      {
        symbol_kind_type yykind = yysym.kind ();
        yyo << (yykind < YYNTOKENS ? "token" : "nterm")
            << ' ' << yysym.name () << " ("
            << yysym.location << ": ";
        switch (yykind)
    {
      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ExpressionOperator > (); }
#line 567 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ExpressionOperator > (); }
#line 573 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ExpressionOperator > (); }
#line 579 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ExpressionOperator > (); }
#line 585 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ExpressionOperator > (); }
#line 591 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 597 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 603 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 609 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 615 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 621 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 627 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 633 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 639 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 645 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 651 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 657 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 663 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 669 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 675 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 681 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 687 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 693 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::FilterOperation > (); }
#line 699 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < std::string > (); }
#line 705 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_INFO_FIELD: // "info key"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < std::string > (); }
#line 711 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < int > (); }
#line 717 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < double > (); }
#line 723 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_BOOL_VALUE: // "boolean value"
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < bool > (); }
#line 729 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_filter: // filter
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < bool > (); }
#line 735 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_simple_expr: // simple_expr
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < bool > (); }
#line 741 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_any_op: // any_op
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::ComparisonOperator > (); }
#line 747 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_str_op: // str_op
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::StringOperator > (); }
#line 753 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_num_op: // num_op
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::NumericalOperator > (); }
#line 759 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_bool_op: // bool_op
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < kim::VariantFilter::BooleanOperator > (); }
#line 765 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_string: // string
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < std::string > (); }
#line 771 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_integer: // integer
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < int > (); }
#line 777 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_number: // number
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < double > (); }
#line 783 "variant_filter_parser.cpp"
        break;

      case symbol_kind::S_boolean: // boolean
#line 91 "variant_filter_parser.ypp"
                 { yyo << yysym.value.template as < bool > (); }
#line 789 "variant_filter_parser.cpp"
        break;

      default:
        break;
    }
        yyo << ')';
      }
  }
#endif

  void
  VariantFilterParser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  VariantFilterParser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  VariantFilterParser::yypop_ (int n) YY_NOEXCEPT
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  VariantFilterParser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  VariantFilterParser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  VariantFilterParser::debug_level_type
  VariantFilterParser::debug_level () const
  {
    return yydebug_;
  }

  void
  VariantFilterParser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  VariantFilterParser::state_type
  VariantFilterParser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool
  VariantFilterParser::yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  VariantFilterParser::yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int
  VariantFilterParser::operator() ()
  {
    return parse ();
  }

  int
  VariantFilterParser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // Discard the LAC context in case there still is one left from a
    // previous invocation.
    yy_lac_discard_ ("init");

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';
    YY_STACK_PRINT ();

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            symbol_type yylookahead (yylex (driver));
            yyla.move (yylookahead);
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    if (yyla.kind () == symbol_kind::S_YYerror)
    {
      // The scanner already issued an error message, process directly
      // to error recovery.  But do not keep the error token as
      // lookahead, it is too special and may lead us to an endless
      // loop in error recovery. */
      yyla.kind_ = symbol_kind::S_YYUNDEF;
      goto yyerrlab1;
    }

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.kind ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind ())
      {
        if (!yy_lac_establish_ (yyla.kind ()))
          goto yyerrlab;
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        if (!yy_lac_establish_ (yyla.kind ()))
          goto yyerrlab;

        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    yy_lac_discard_ ("shift");
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
      switch (yyr1_[yyn])
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        yylhs.value.emplace< bool > ();
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        yylhs.value.emplace< double > ();
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        yylhs.value.emplace< int > ();
        break;

      case symbol_kind::S_bool_op: // bool_op
        yylhs.value.emplace< kim::VariantFilter::BooleanOperator > ();
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        yylhs.value.emplace< kim::VariantFilter::ComparisonOperator > ();
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        yylhs.value.emplace< kim::VariantFilter::ExpressionOperator > ();
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        yylhs.value.emplace< kim::VariantFilter::FilterOperation > ();
        break;

      case symbol_kind::S_num_op: // num_op
        yylhs.value.emplace< kim::VariantFilter::NumericalOperator > ();
        break;

      case symbol_kind::S_str_op: // str_op
        yylhs.value.emplace< kim::VariantFilter::StringOperator > ();
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        yylhs.value.emplace< std::string > ();
        break;

      default:
        break;
    }


      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 2: // expr: filter
#line 97 "variant_filter_parser.ypp"
             { driver._result = yystack_[0].value.as < bool > (); }
#line 1120 "variant_filter_parser.cpp"
    break;

  case 3: // filter: filter "conjunction operator ('AND' or '&&')" filter
#line 100 "variant_filter_parser.ypp"
                      { yylhs.value.as < bool > () = yystack_[2].value.as < bool > () && yystack_[0].value.as < bool > (); }
#line 1126 "variant_filter_parser.cpp"
    break;

  case 4: // filter: filter "disjunction operator ('or' or '||')" filter
#line 101 "variant_filter_parser.ypp"
                     { yylhs.value.as < bool > () = yystack_[2].value.as < bool > () || yystack_[0].value.as < bool > (); }
#line 1132 "variant_filter_parser.cpp"
    break;

  case 5: // filter: "negation operator ('not' or '!')" filter
#line 102 "variant_filter_parser.ypp"
               { yylhs.value.as < bool > () = !yystack_[0].value.as < bool > (); }
#line 1138 "variant_filter_parser.cpp"
    break;

  case 6: // filter: "opening parenthesis '('" filter "closing parenthesis ')'"
#line 103 "variant_filter_parser.ypp"
                                   { yylhs.value.as < bool > () = yystack_[1].value.as < bool > (); }
#line 1144 "variant_filter_parser.cpp"
    break;

  case 7: // filter: simple_expr
#line 104 "variant_filter_parser.ypp"
                { yylhs.value.as < bool > () = yystack_[0].value.as < bool > (); }
#line 1150 "variant_filter_parser.cpp"
    break;

  case 8: // simple_expr: "field attribute 'CHROM'" str_op string
#line 108 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onChrom(yystack_[1].value.as < kim::VariantFilter::StringOperator > (), yystack_[0].value.as < std::string > ()); }
#line 1156 "variant_filter_parser.cpp"
    break;

  case 9: // simple_expr: "field attribute 'POS'" num_op integer
#line 109 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onPos(yystack_[1].value.as < kim::VariantFilter::NumericalOperator > (), yystack_[0].value.as < int > ()); }
#line 1162 "variant_filter_parser.cpp"
    break;

  case 10: // simple_expr: "field attribute 'ID'" str_op string
#line 110 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onID(yystack_[1].value.as < kim::VariantFilter::StringOperator > (), yystack_[0].value.as < std::string > ()); }
#line 1168 "variant_filter_parser.cpp"
    break;

  case 11: // simple_expr: "field attribute 'QUAL'" num_op number
#line 111 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onQuality(yystack_[1].value.as < kim::VariantFilter::NumericalOperator > (), yystack_[0].value.as < double > ()); }
#line 1174 "variant_filter_parser.cpp"
    break;

  case 12: // simple_expr: "field attribute 'FILTER'" str_op string
#line 112 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onFilter(yystack_[1].value.as < kim::VariantFilter::StringOperator > (), yystack_[0].value.as < std::string > ()); }
#line 1180 "variant_filter_parser.cpp"
    break;

  case 13: // simple_expr: "field attribute 'INFO'" "info key" any_op string
#line 113 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onInfo(yystack_[2].value.as < std::string > (), yystack_[1].value.as < kim::VariantFilter::ComparisonOperator > (), yystack_[0].value.as < std::string > ()); }
#line 1186 "variant_filter_parser.cpp"
    break;

  case 14: // simple_expr: "status 'Ploidy'" num_op integer
#line 114 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onPloidy(yystack_[1].value.as < kim::VariantFilter::NumericalOperator > (), yystack_[0].value.as < int > ()); }
#line 1192 "variant_filter_parser.cpp"
    break;

  case 15: // simple_expr: "status 'SNP' (Single Nucleotide Polymorphism)" bool_op boolean
#line 115 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onSNP(yystack_[1].value.as < kim::VariantFilter::BooleanOperator > (), yystack_[0].value.as < bool > ()); }
#line 1198 "variant_filter_parser.cpp"
    break;

  case 16: // simple_expr: "status 'MSNP' (multi-allelic SNP)" bool_op boolean
#line 116 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onMultiAllelicSNP(yystack_[1].value.as < kim::VariantFilter::BooleanOperator > (), yystack_[0].value.as < bool > ()); }
#line 1204 "variant_filter_parser.cpp"
    break;

  case 17: // simple_expr: "status 'SV' (Structural Variant)" bool_op boolean
#line 117 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onSV(yystack_[1].value.as < kim::VariantFilter::BooleanOperator > (), yystack_[0].value.as < bool > ()); }
#line 1210 "variant_filter_parser.cpp"
    break;

  case 18: // simple_expr: "status 'InDel' (Insertion-Deletion)" bool_op boolean
#line 118 "variant_filter_parser.ypp"
                                            { yylhs.value.as < bool > () = driver.onIndel(yystack_[1].value.as < kim::VariantFilter::BooleanOperator > (), yystack_[0].value.as < bool > ()); }
#line 1216 "variant_filter_parser.cpp"
    break;

  case 19: // any_op: "equality operator ('=')"
#line 122 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1222 "variant_filter_parser.cpp"
    break;

  case 20: // any_op: "regular expression match operator ('~')"
#line 123 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1228 "variant_filter_parser.cpp"
    break;

  case 21: // any_op: "inequality operator ('!=' or '<>')"
#line 124 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1234 "variant_filter_parser.cpp"
    break;

  case 22: // any_op: "order relationship operator '<'"
#line 125 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1240 "variant_filter_parser.cpp"
    break;

  case 23: // any_op: "order relationship operator '<='"
#line 126 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1246 "variant_filter_parser.cpp"
    break;

  case 24: // any_op: "order relationship operator '>'"
#line 127 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1252 "variant_filter_parser.cpp"
    break;

  case 25: // any_op: "order relationship operator '>='"
#line 128 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::ComparisonOperator > () = yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > (); }
#line 1258 "variant_filter_parser.cpp"
    break;

  case 26: // str_op: "equality operator ('=')"
#line 132 "variant_filter_parser.ypp"
               { yylhs.value.as < kim::VariantFilter::StringOperator > () = VariantFilter::StringOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1264 "variant_filter_parser.cpp"
    break;

  case 27: // str_op: "regular expression match operator ('~')"
#line 133 "variant_filter_parser.ypp"
               { yylhs.value.as < kim::VariantFilter::StringOperator > () = VariantFilter::StringOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1270 "variant_filter_parser.cpp"
    break;

  case 28: // str_op: "inequality operator ('!=' or '<>')"
#line 134 "variant_filter_parser.ypp"
               { yylhs.value.as < kim::VariantFilter::StringOperator > () = VariantFilter::StringOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1276 "variant_filter_parser.cpp"
    break;

  case 29: // num_op: "equality operator ('=')"
#line 138 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1282 "variant_filter_parser.cpp"
    break;

  case 30: // num_op: "inequality operator ('!=' or '<>')"
#line 139 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1288 "variant_filter_parser.cpp"
    break;

  case 31: // num_op: "order relationship operator '<'"
#line 140 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1294 "variant_filter_parser.cpp"
    break;

  case 32: // num_op: "order relationship operator '<='"
#line 141 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1300 "variant_filter_parser.cpp"
    break;

  case 33: // num_op: "order relationship operator '>'"
#line 142 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1306 "variant_filter_parser.cpp"
    break;

  case 34: // num_op: "order relationship operator '>='"
#line 143 "variant_filter_parser.ypp"
                      { yylhs.value.as < kim::VariantFilter::NumericalOperator > () = VariantFilter::NumericalOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1312 "variant_filter_parser.cpp"
    break;

  case 35: // bool_op: "equality operator ('=')"
#line 147 "variant_filter_parser.ypp"
               { yylhs.value.as < kim::VariantFilter::BooleanOperator > () = VariantFilter::BooleanOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1318 "variant_filter_parser.cpp"
    break;

  case 36: // bool_op: "inequality operator ('!=' or '<>')"
#line 148 "variant_filter_parser.ypp"
               { yylhs.value.as < kim::VariantFilter::BooleanOperator > () = VariantFilter::BooleanOperator(yystack_[0].value.as < kim::VariantFilter::ComparisonOperator > ()); }
#line 1324 "variant_filter_parser.cpp"
    break;

  case 37: // string: "string value"
#line 152 "variant_filter_parser.ypp"
                { yylhs.value.as < std::string > () = yystack_[0].value.as < std::string > (); }
#line 1330 "variant_filter_parser.cpp"
    break;

  case 38: // string: "integer value"
#line 153 "variant_filter_parser.ypp"
                { yylhs.value.as < std::string > () = std::to_string(yystack_[0].value.as < int > ()); }
#line 1336 "variant_filter_parser.cpp"
    break;

  case 39: // string: "number value"
#line 154 "variant_filter_parser.ypp"
                { yylhs.value.as < std::string > () = std::to_string(yystack_[0].value.as < double > ()); }
#line 1342 "variant_filter_parser.cpp"
    break;

  case 40: // string: "boolean value"
#line 155 "variant_filter_parser.ypp"
                { yylhs.value.as < std::string > () = yystack_[0].value.as < bool > () ? "true" : "false"; }
#line 1348 "variant_filter_parser.cpp"
    break;

  case 41: // integer: "integer value"
#line 158 "variant_filter_parser.ypp"
                    { yylhs.value.as < int > () = yystack_[0].value.as < int > (); }
#line 1354 "variant_filter_parser.cpp"
    break;

  case 42: // number: "integer value"
#line 161 "variant_filter_parser.ypp"
               { yylhs.value.as < double > () = yystack_[0].value.as < int > (); }
#line 1360 "variant_filter_parser.cpp"
    break;

  case 43: // number: "number value"
#line 162 "variant_filter_parser.ypp"
               { yylhs.value.as < double > () = yystack_[0].value.as < double > (); }
#line 1366 "variant_filter_parser.cpp"
    break;

  case 44: // boolean: "boolean value"
#line 166 "variant_filter_parser.ypp"
               { yylhs.value.as < bool > () = yystack_[0].value.as < bool > (); }
#line 1372 "variant_filter_parser.cpp"
    break;

  case 45: // boolean: "integer value"
#line 167 "variant_filter_parser.ypp"
               { yylhs.value.as < bool > () = yystack_[0].value.as < int > (); }
#line 1378 "variant_filter_parser.cpp"
    break;


#line 1382 "variant_filter_parser.cpp"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        context yyctx (*this, yyla);
        std::string msg = yysyntax_error_ (yyctx);
        error (yyla.location, YY_MOVE (msg));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.kind () == symbol_kind::S_YYEOF)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    // Pop stack until we find a state that shifts the error token.
    for (;;)
      {
        yyn = yypact_[+yystack_[0].state];
        if (!yy_pact_value_is_default_ (yyn))
          {
            yyn += symbol_kind::S_YYerror;
            if (0 <= yyn && yyn <= yylast_
                && yycheck_[yyn] == symbol_kind::S_YYerror)
              {
                yyn = yytable_[yyn];
                if (0 < yyn)
                  break;
              }
          }

        // Pop the current state because it cannot handle the error token.
        if (yystack_.size () == 1)
          YYABORT;

        yyerror_range[1].location = yystack_[0].location;
        yy_destroy_ ("Error: popping", yystack_[0]);
        yypop_ ();
        YY_STACK_PRINT ();
      }
    {
      stack_symbol_type error_token;

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      yy_lac_discard_ ("error recovery");
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
    goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
  `-------------------------------------*/
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;


  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
  `-----------------------------------*/
  yyabortlab:
    yyresult = 1;
    goto yyreturn;


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    YY_STACK_PRINT ();
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  VariantFilterParser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

  const char *
  VariantFilterParser::symbol_name (symbol_kind_type yysymbol)
  {
    static const char *const yy_sname[] =
    {
    "end of file", "error", "invalid token",
  "conjunction operator ('AND' or '&&')",
  "disjunction operator ('or' or '||')",
  "negation operator ('not' or '!')", "opening parenthesis '('",
  "closing parenthesis ')'", "equality operator ('=')",
  "regular expression match operator ('~')",
  "inequality operator ('!=' or '<>')", "order relationship operator '<'",
  "order relationship operator '<='", "order relationship operator '>'",
  "order relationship operator '>='", "field attribute 'CHROM'",
  "field attribute 'POS'", "field attribute 'ID'",
  "field attribute 'QUAL'", "field attribute 'FILTER'",
  "field attribute 'INFO'", "status 'Ploidy'",
  "status 'SNP' (Single Nucleotide Polymorphism)",
  "status 'MSNP' (multi-allelic SNP)", "status 'SV' (Structural Variant)",
  "status 'InDel' (Insertion-Deletion)", "string value", "info key",
  "integer value", "number value", "boolean value", "$accept", "expr",
  "filter", "simple_expr", "any_op", "str_op", "num_op", "bool_op",
  "string", "integer", "number", "boolean", YY_NULLPTR
    };
    return yy_sname[yysymbol];
  }



  // VariantFilterParser::context.
  VariantFilterParser::context::context (const VariantFilterParser& yyparser, const symbol_type& yyla)
    : yyparser_ (yyparser)
    , yyla_ (yyla)
  {}

  int
  VariantFilterParser::context::expected_tokens (symbol_kind_type yyarg[], int yyargn) const
  {
    // Actual number of expected tokens
    int yycount = 0;

#if YYDEBUG
    // Execute LAC once. We don't care if it is successful, we
    // only do it for the sake of debugging output.
    if (!yyparser_.yy_lac_established_)
      yyparser_.yy_lac_check_ (yyla_.kind ());
#endif

    for (int yyx = 0; yyx < YYNTOKENS; ++yyx)
      {
        symbol_kind_type yysym = YY_CAST (symbol_kind_type, yyx);
        if (yysym != symbol_kind::S_YYerror
            && yysym != symbol_kind::S_YYUNDEF
            && yyparser_.yy_lac_check_ (yysym))
          {
            if (!yyarg)
              ++yycount;
            else if (yycount == yyargn)
              return 0;
            else
              yyarg[yycount++] = yysym;
          }
      }
    if (yyarg && yycount == 0 && 0 < yyargn)
      yyarg[0] = symbol_kind::S_YYEMPTY;
    return yycount;
  }




  bool
  VariantFilterParser::yy_lac_check_ (symbol_kind_type yytoken) const
  {
    // Logically, the yylac_stack's lifetime is confined to this function.
    // Clear it, to get rid of potential left-overs from previous call.
    yylac_stack_.clear ();
    // Reduce until we encounter a shift and thereby accept the token.
#if YYDEBUG
    YYCDEBUG << "LAC: checking lookahead " << symbol_name (yytoken) << ':';
#endif
    std::ptrdiff_t lac_top = 0;
    while (true)
      {
        state_type top_state = (yylac_stack_.empty ()
                                ? yystack_[lac_top].state
                                : yylac_stack_.back ());
        int yyrule = yypact_[+top_state];
        if (yy_pact_value_is_default_ (yyrule)
            || (yyrule += yytoken) < 0 || yylast_ < yyrule
            || yycheck_[yyrule] != yytoken)
          {
            // Use the default action.
            yyrule = yydefact_[+top_state];
            if (yyrule == 0)
              {
                YYCDEBUG << " Err\n";
                return false;
              }
          }
        else
          {
            // Use the action from yytable.
            yyrule = yytable_[yyrule];
            if (yy_table_value_is_error_ (yyrule))
              {
                YYCDEBUG << " Err\n";
                return false;
              }
            if (0 < yyrule)
              {
                YYCDEBUG << " S" << yyrule << '\n';
                return true;
              }
            yyrule = -yyrule;
          }
        // By now we know we have to simulate a reduce.
        YYCDEBUG << " R" << yyrule - 1;
        // Pop the corresponding number of values from the stack.
        {
          std::ptrdiff_t yylen = yyr2_[yyrule];
          // First pop from the LAC stack as many tokens as possible.
          std::ptrdiff_t lac_size = std::ptrdiff_t (yylac_stack_.size ());
          if (yylen < lac_size)
            {
              yylac_stack_.resize (std::size_t (lac_size - yylen));
              yylen = 0;
            }
          else if (lac_size)
            {
              yylac_stack_.clear ();
              yylen -= lac_size;
            }
          // Only afterwards look at the main stack.
          // We simulate popping elements by incrementing lac_top.
          lac_top += yylen;
        }
        // Keep top_state in sync with the updated stack.
        top_state = (yylac_stack_.empty ()
                     ? yystack_[lac_top].state
                     : yylac_stack_.back ());
        // Push the resulting state of the reduction.
        state_type state = yy_lr_goto_state_ (top_state, yyr1_[yyrule]);
        YYCDEBUG << " G" << int (state);
        yylac_stack_.push_back (state);
      }
  }

  // Establish the initial context if no initial context currently exists.
  bool
  VariantFilterParser::yy_lac_establish_ (symbol_kind_type yytoken)
  {
    /* Establish the initial context for the current lookahead if no initial
       context is currently established.

       We define a context as a snapshot of the parser stacks.  We define
       the initial context for a lookahead as the context in which the
       parser initially examines that lookahead in order to select a
       syntactic action.  Thus, if the lookahead eventually proves
       syntactically unacceptable (possibly in a later context reached via a
       series of reductions), the initial context can be used to determine
       the exact set of tokens that would be syntactically acceptable in the
       lookahead's place.  Moreover, it is the context after which any
       further semantic actions would be erroneous because they would be
       determined by a syntactically unacceptable token.

       yy_lac_establish_ should be invoked when a reduction is about to be
       performed in an inconsistent state (which, for the purposes of LAC,
       includes consistent states that don't know they're consistent because
       their default reductions have been disabled).

       For parse.lac=full, the implementation of yy_lac_establish_ is as
       follows.  If no initial context is currently established for the
       current lookahead, then check if that lookahead can eventually be
       shifted if syntactic actions continue from the current context.  */
    if (yy_lac_established_)
      return true;
    else
      {
#if YYDEBUG
        YYCDEBUG << "LAC: initial context established for "
                 << symbol_name (yytoken) << '\n';
#endif
        yy_lac_established_ = true;
        return yy_lac_check_ (yytoken);
      }
  }

  // Discard any previous initial lookahead context.
  void
  VariantFilterParser::yy_lac_discard_ (const char* event)
  {
   /* Discard any previous initial lookahead context because of Event,
      which may be a lookahead change or an invalidation of the currently
      established initial context for the current lookahead.

      The most common example of a lookahead change is a shift.  An example
      of both cases is syntax error recovery.  That is, a syntax error
      occurs when the lookahead is syntactically erroneous for the
      currently established initial context, so error recovery manipulates
      the parser stacks to try to find a new initial context in which the
      current lookahead is syntactically acceptable.  If it fails to find
      such a context, it discards the lookahead.  */
    if (yy_lac_established_)
      {
        YYCDEBUG << "LAC: initial context discarded due to "
                 << event << '\n';
        yy_lac_established_ = false;
      }
  }


  int
  VariantFilterParser::yy_syntax_error_arguments_ (const context& yyctx,
                                                 symbol_kind_type yyarg[], int yyargn) const
  {
    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
         In the first two cases, it might appear that the current syntax
         error should have been detected in the previous state when
         yy_lac_check was invoked.  However, at that time, there might
         have been a different syntax error that discarded a different
         initial context during error recovery, leaving behind the
         current lookahead.
    */

    if (!yyctx.lookahead ().empty ())
      {
        if (yyarg)
          yyarg[0] = yyctx.token ();
        int yyn = yyctx.expected_tokens (yyarg ? yyarg + 1 : yyarg, yyargn - 1);
        return yyn + 1;
      }
    return 0;
  }

  // Generate an error message.
  std::string
  VariantFilterParser::yysyntax_error_ (const context& yyctx) const
  {
    // Its maximum.
    enum { YYARGS_MAX = 5 };
    // Arguments of yyformat.
    symbol_kind_type yyarg[YYARGS_MAX];
    int yycount = yy_syntax_error_arguments_ (yyctx, yyarg, YYARGS_MAX);

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
      default: // Avoid compiler warnings.
        YYCASE_ (0, YY_("syntax error"));
        YYCASE_ (1, YY_("syntax error, unexpected %s"));
        YYCASE_ (2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_ (3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_ (4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_ (5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    std::ptrdiff_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += symbol_name (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char VariantFilterParser::yypact_ninf_ = -29;

  const signed char VariantFilterParser::yytable_ninf_ = -1;

  const signed char
  VariantFilterParser::yypact_[] =
  {
       0,     0,     0,     4,    35,     4,    35,     4,   -24,    35,
      48,    48,    48,    48,     8,    36,   -29,   -29,    30,   -29,
     -29,   -29,   -19,   -29,   -29,   -29,   -29,   -29,   -29,    16,
     -19,    37,   -19,    18,    16,   -29,   -29,    31,    31,    31,
      31,   -29,     0,     0,   -29,   -29,   -29,   -29,   -29,   -29,
     -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,
     -29,   -29,   -29,   -29,   -19,   -29,   -29,   -29,   -29,   -29,
     -29,   -29,   -29,    54,   -29
  };

  const signed char
  VariantFilterParser::yydefact_[] =
  {
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     7,     5,     0,    26,
      27,    28,     0,    29,    30,    31,    32,    33,    34,     0,
       0,     0,     0,     0,     0,    35,    36,     0,     0,     0,
       0,     1,     0,     0,     6,    37,    38,    39,    40,     8,
      41,     9,    10,    42,    43,    11,    12,    19,    20,    21,
      22,    23,    24,    25,     0,    14,    45,    44,    15,    16,
      17,    18,     3,     4,    13
  };

  const signed char
  VariantFilterParser::yypgoto_[] =
  {
     -29,   -29,    -1,   -29,   -29,    55,    29,    39,   -28,    33,
     -29,    15
  };

  const signed char
  VariantFilterParser::yydefgoto_[] =
  {
       0,    14,    15,    16,    64,    22,    29,    37,    49,    51,
      55,    68
  };

  const signed char
  VariantFilterParser::yytable_[] =
  {
      17,    18,    52,    33,    56,     1,     2,    45,    41,    46,
      47,    48,    19,    20,    21,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    57,    58,    59,    60,
      61,    62,    63,    42,    43,    31,    74,    44,    34,    42,
      43,    72,    73,    23,    50,    24,    25,    26,    27,    28,
      38,    39,    40,    69,    70,    71,    35,    42,    36,    66,
      30,    67,    32,     0,     0,    53,    54,    65
  };

  const signed char
  VariantFilterParser::yycheck_[] =
  {
       1,     2,    30,    27,    32,     5,     6,    26,     0,    28,
      29,    30,     8,     9,    10,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,     8,     9,    10,    11,
      12,    13,    14,     3,     4,     6,    64,     7,     9,     3,
       4,    42,    43,     8,    28,    10,    11,    12,    13,    14,
      11,    12,    13,    38,    39,    40,     8,     3,    10,    28,
       5,    30,     7,    -1,    -1,    28,    29,    34
  };

  const signed char
  VariantFilterParser::yystos_[] =
  {
       0,     5,     6,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    32,    33,    34,    33,    33,     8,
       9,    10,    36,     8,    10,    11,    12,    13,    14,    37,
      36,    37,    36,    27,    37,     8,    10,    38,    38,    38,
      38,     0,     3,     4,     7,    26,    28,    29,    30,    39,
      28,    40,    39,    28,    29,    41,    39,     8,     9,    10,
      11,    12,    13,    14,    35,    40,    28,    30,    42,    42,
      42,    42,    33,    33,    39
  };

  const signed char
  VariantFilterParser::yyr1_[] =
  {
       0,    31,    32,    33,    33,    33,    33,    33,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    35,
      35,    35,    35,    35,    35,    35,    36,    36,    36,    37,
      37,    37,    37,    37,    37,    38,    38,    39,    39,    39,
      39,    40,    41,    41,    42,    42
  };

  const signed char
  VariantFilterParser::yyr2_[] =
  {
       0,     2,     1,     3,     3,     2,     3,     1,     3,     3,
       3,     3,     3,     4,     3,     3,     3,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1
  };




#if YYDEBUG
  const unsigned char
  VariantFilterParser::yyrline_[] =
  {
       0,    97,    97,   100,   101,   102,   103,   104,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   122,
     123,   124,   125,   126,   127,   128,   132,   133,   134,   138,
     139,   140,   141,   142,   143,   147,   148,   152,   153,   154,
     155,   158,   161,   162,   166,   167
  };

  void
  VariantFilterParser::yy_stack_print_ () const
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  void
  VariantFilterParser::yy_reduce_print_ (int yyrule) const
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG


#line 15 "variant_filter_parser.ypp"
} // kim
#line 1997 "variant_filter_parser.cpp"

#line 169 "variant_filter_parser.ypp"


void kim::VariantFilterParser::error(const location_type &l, const std::string &msg) {
  driver._error_message = msg;
  driver.scanner_location = l;
}
