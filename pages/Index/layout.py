import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table

PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"

layout = html.Div([
    dcc.Location(id='url', refresh=False),
    # Use a row instead of a number and add cute buttons
    dbc.Navbar([
    #dbc.Container([
        html.Img(src=PLOTLY_LOGO, height="40px"),
        html.P("...", style={"color": "transparent"}),
        #dbc.NavItem(dbc.NavbarBrand("PlantApp")),
        #dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),

        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("About", id="home_about", external_link=True, href="/about"),
            dbc.DropdownMenuItem("Login", id="home_login", external_link=True, href="/login"),
            ], label="Home", color="primary", className="m-1",
        ),
        
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Available DBs", id="available_dbs", external_link=True, href="/SQNce/available_dbs"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Descriptions", header=True),
            dbc.DropdownMenuItem("Annotations", id="annotations", external_link=True, href="/SQNce/annotations"),
            dbc.DropdownMenuItem("Gene Symbols", id="symbols", disabled=True),
            dbc.DropdownMenuItem("Get Genes from Family Names", id="families_familyIDs", external_link=True, href="/SQNce/families_familyIDs"),
            dbc.DropdownMenuItem("Get Family Names of Genes", id="families_geneIDs", external_link=True, href="/SQNce/families_geneIDs"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Relations", header=True),
            dbc.DropdownMenuItem("Blast Best Hits (BBHs)", id="BBHs", external_link=True, href="/SQNce/BBHs"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Coordinates", header=True),
            dbc.DropdownMenuItem("Coordinates", id="coordinates", external_link=True, href="/SQNce/coordinates"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Sequences", header=True),
            dbc.DropdownMenuItem("Sequences", id="proteins", external_link=True, href="/SQNce/proteins"),
            dbc.DropdownMenuItem("Promoters", id="promoters", external_link=True, href="/SQNce/promoters"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Transcriptomic Data", header=True),
            dbc.DropdownMenuItem("Transcriptomic Meta-data", header=True),
            dbc.DropdownMenuItem("Omics", id="omics", external_link=True, href="/SQNce/omics"),
            ], label="SQNce", color="secondary", className="m-1",
            toggle_style={"background": "transparent","border-color": "#f8f9fa"},
        ),
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Simple Tree", id="simple_tree", external_link=True, href="/apps/simple_tree"),
            dbc.DropdownMenuItem("Genome Graph", id="genome_graph", external_link=True, href="/apps/genome_graph"),
            dbc.DropdownMenuItem("GO Enrichment", id="go_enrichment", disabled=True),
            dbc.DropdownMenuItem("Mapping Summary", id="mapping_summary", external_link=True, href="/apps/mapping_summary"),
            ], label="Apps", color="secondary", className="m-1",
            toggle_style={"background": "transparent","border-color": "#f8f9fa"},
        ),
        dbc.Button("Downloads", id="downloads", outline=True, color="secondary",
                    className="m-1", disabled=True),
        dbc.Button("Feedback", id="feedback", outline=True, color="secondary",
                    className="m-1", disabled=True),
        ],
        #align="start",
        #style={"margin-left": "15px"}
        dark =True,
        color="dark",
        sticky = True,
    ),

    dbc.Row(
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
])