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

        dbc.Button("About", id="about", outline=False, color="primary", className="m-1", external_link=True, href="/about"),
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
            ], label="SQNce", color="primary", className="m-1",
            toggle_style={"background": "transparent","border-color": "#f8f9fa"},
        ),
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Simple Tree", id="simple_tree", external_link=True, href="/apps/simple_tree"),
            dbc.DropdownMenuItem("GO Enrichment", id="go_enrichment", disabled=True),
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
        dbc.Col(
            html.Div(
                #html.H1("Welcome to PlantApp"),
                #html.P('Dash converts Python classes into HTML'),
                dcc.Markdown('''
                    ## Welcome to PlantApp

                    PlantApp contains different tools for comparative genomics analysis and it uses [SQNce](https://github.com/eporetsky/SQNce),
                    as a back-end. SQNce is a SQLite-based database framework for parsing commonly used gene-annotation files.
                    PlantApp aggregates publicly available data from different sequenced plant species and genotypes and offers an
                    easy-to-access and easy-to-download database.
                    ''')
                ),
            width=6),

        dbc.Col(
            html.Div(
                dcc.Markdown('''
                    ## Available data

                    Will be added later.
                '''
                ),
            ),
        width=6),
        ]
    ),
    #html.Div(id="tab_content"),


      #dbc.Row([table])
    ],
)