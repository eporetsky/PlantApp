import sqlite3
import os 
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table




def register_callbacks(dashapp):
    from dash.dependencies import Input, Output

    import pandas as pd
    from collections import OrderedDict
    import dash_bootstrap_components as dbc
    from dash import Dash, dcc, html, Input, Output, State, dash_table

    import csv

    from flask import Flask
    import sqlite3
    import numpy as np
    import plotly.graph_objs as go

    import os
    import zlib


    # Set the cwd and SQNce.db path to be platform independent  
    if "plantapp" not in os.getcwd().lower():
        cwd = "/home/eporetsky/plantapp" # for server hosting
    else:
        cwd = os.getcwd() # for personal computer
    sqnce_path = os.path.join(cwd, "SQNce.db")
    
    @dashapp.callback(Output('tab_content', 'children'),
                      Input('url', 'pathname'))
    def display_page_content(pathname):
        #print(session.get('username', None))
        print(pathname)
        if pathname == "/index/" or pathname == "/":
            return tab_about

        pathname = pathname.split("/")[-1]
        print(pathname)
        if pathname == "index":
            return tab_about
        if pathname == "changelog":
            return tab_changelog

    tab_about = dbc.Row(
        [

        dbc.Col([
            dbc.Container([
                html.H1("Welcome to PlantApp", className="display-6"),
                html.P("""PlantApp contains different tools for comparative genomics analysis and it uses SQNce,
                    as a back-end. SQNce is a SQLite-based database framework for parsing commonly used gene-annotation files.
                    PlantApp aggregates publicly available data from different sequenced plant species and genotypes and offers an
                    easy-to-access and easy-to-download database.""",
                    #"featured content or information.",
                    className="lead",
                ),
                html.Hr(className="my-2"),
                html.P("""The website, database queries and apps are all under development but all the individual components function are functional.  
                """
                ),
                ],
                fluid=True,
                className="py-3"
                ),
            ], width={"size": 4, "order":1, "offset": 1},
        ),


        dbc.Col([
            dbc.Container([
                html.H1("SQNce Queries", className="display-6"),
                html.P("""PlantApp uses built-in queries to quickly get user-requested data. The expandable menu below provides a quick
                    overview of the available queries and how to use them. Queries are added based on existing needs and will be periodically updated.
                    """,
                    className="lead",
                ),
                html.Hr(className="my-2"),
                
                dbc.Accordion([
                    dbc.AccordionItem(
                        [
                            html.P("""PlantApp uses a SQNce database file to store common types of gene-related information from multiple species.
                                Different information from different species will be continuously added. Press on the Available Data button to see
                                which genomes where included for each type of information, including the source of the parsed data."""),
                        ],
                        title="Available Data",
                    ),
                    dbc.AccordionItem(
                        [
                            html.P("""Multiple types of gene annotations provide useful information about known and predicted gene functions.
                                The purpose of this tab is to integrate multiple sources of information for a large number of genomes in a single
                                convenient place. 
                                """),
                        ],
                        title="Gene Description",
                    ),
                    dbc.AccordionItem(
                        [
                            html.P("""Comparative genetic approaches include the understanding of relations between genes. This section is mainly 
                                focused on sequence similarity as a convenient way to identify gene families and homologs within and across species.
                                Sequence similarity can be summarized and analyzed in multiple way, some of which are included here.
                                """),
                        ],
                        title="Gene Relations",
                    ),
                ], start_collapsed=True,)


                ],
                fluid=True,
                className="py-3"
            ),

            dbc.Container([
                html.H1("Apps", className="display-6"),
                html.P("""PlantApp includes a number of simple apps that could take as an input some of the query results and provide additional
                    information.
                    """,
                    className="lead",
                ),
                html.Hr(className="my-2"),
                
                dbc.Accordion([
                    dbc.AccordionItem(
                        [
                            html.P("""A simple tree builder that takes gene IDs and aligns their primary protein sequences and generates a phylogeny tree."""),
                        ],
                        title="Simple Tree",
                    ),
                    dbc.AccordionItem(
                        [
                            html.P("""In progress.
                                """),
                        ],
                        title="Protein Domains",
                    ),
                    dbc.AccordionItem(
                        [
                            html.P("""In progress.
                                """),
                        ],
                        title="Gene Ontology Enrichment",
                    ),
                ], start_collapsed=True,)


                ],
                fluid=True,
                className="py-3"
            ),
            ], width={"size": 4, "order":2, "offset": 2},
        ),
    ])

    tab_changelog = dbc.Row([
        dbc.Col([
        dbc.Container([
            html.H1("Changelog", className="display-6"),
            html.P("""The changelog section contains notes on when different sections were added to PlantApp.
                      In the future it is likely to include the different datasets that were added to the SQNce database.""",
                #"featured content or information.",
                className="lead",
            ),
            html.Hr(className="my-2"),

        dbc.Accordion(
        [
            dbc.AccordionItem(
                "1. Added the changelog section", title="March 2022"
            ),
            dbc.AccordionItem(
                [html.P("1. Added GO enrichment analysis"),
                 html.P("2. Added option to identify mapping candidates")], 
                 title="May 2022"),
            dbc.AccordionItem(
                [html.P("1. Added the gene symbols tab"),
                 html.P("2. Added a primitive login page")], 
                 title="March 2022"),
            dbc.AccordionItem(
                [html.P("1. Added genome grapth from gene list app"),
                 html.P("2. Improved index page with dash-bootstrap-components")], 
                 title="February 2022"),
            dbc.AccordionItem(
                [html.P("1. Added genome grapth from gene list app"),
                 html.P("2. Improved index page with dash-bootstrap-components")], 
                 title="January 2022"),
            dbc.AccordionItem(
                [html.P("1. Added the simpleTree app"),
                 html.P("2. Added family annotations best on interpro"),
                 html.P("3. Added gene candidate annotation using genomic coordinates"),
                 html.P("4. Added the menu tabs"),
                 html.P("5. Embeded the different Dash apps into one Flask apps")], 
                 title="December 2021"),
            dbc.AccordionItem(
                [html.P("1. Added promoter based sequences")], 
                 title="November 2021"),
            dbc.AccordionItem(
                [html.P("1. Added protein sequences"),
                 html.P("2. Added gene annotations"), 
                 html.P("3. Added query for available datasets"), ], 
                 title="October 2021"),
            dbc.AccordionItem(
                [html.P("1. Initial PlantApp commit")], 
                 title="August 2021"),
        ],
        start_collapsed=False,
        ), 
    ], fluid=True, className="py-3", ),
    ], width={"size": 6, "order":1, "offset": 3})
    ])