import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table

PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"

layout = dbc.Container(
    [
        dcc.Store(id="store"),

        #html.H1("Dynamically rendered tab content"),
        #html.Hr(),

    dbc.Navbar(
        children = [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("PlantApp", className="ml-2")),
                    dbc.Col(dbc.NavLink("About", href="/about/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("SQNce", href="/SQNce/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Apps", href="/apps/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Downloads", href="/downloads/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Feedback", href="/feedback/", className="ml-2", external_link=True)),
                    ],
                    align="center",# no_gutters=True,
                ),
                href="https://www.plantapp.org",
            ),
            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            #dbc.Collapse(search_bar, id="navbar-collapse", navbar=True, is_open=False),
        ],
    color="dark", dark=True,
    ),

    html.A(dbc.Row([
            dbc.Button("Available DBs", id="available_dbs", className="me-1"),
            dbc.Button("Annotations", id="annotations", className="me-1"),
            dbc.DropdownMenu([dbc.DropdownMenuItem("Using Family Names", id="families_familyIDs"),
                              dbc.DropdownMenuItem("Using Gene IDs", id="families_geneIDs")], label="Families"),
            dbc.Button("Coordinates", id="coordinates", className="me-1"),
            dbc.Button("Sequences", id="proteins", className="me-1"),
            dbc.Button("Promoters", id="promoters", className="me-1"),
            dbc.Button("Omics", id="omics", className="me-1"),
    ])),

 
    html.Div(id="tab_content"),

    ],
    fluid=True  
)