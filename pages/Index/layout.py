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
            [dbc.DropdownMenuItem("About", id="home_about", external_link=True, href="/index"),
            dbc.DropdownMenuItem("Changelog", id="home_changelog", external_link=True, href="/index/changelog"),
            dbc.DropdownMenuItem("Login", id="home_login", external_link=True, href="/login"),
            ], label="Home", color="primary", className="m-1",
        ),
        
        dbc.DropdownMenu(
            [dbc.DropdownMenuItem("Available DBs", id="available_dbs", external_link=True, href="/SQNce/available_dbs"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Gene Descriptions", header=True),
            dbc.DropdownMenuItem("Annotations", id="annotations", external_link=True, href="/SQNce/annotations"),
            dbc.DropdownMenuItem("Gene Symbols", id="symbols", external_link=True, href="/SQNce/symbols"),
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
            dbc.DropdownMenuItem("GO Enrichment", id="go_enrichment", external_link=True, href="/apps/go_enrichment"),
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

    html.Div(id="tab_content"),

    
])