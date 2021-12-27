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
                    dbc.Col(dbc.NavLink("Queries", href="/queries/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Tools", href="/tools/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Downloads", href="/downloads/", className="ml-2", external_link=True)),
                    dbc.Col(dbc.NavLink("Feedback", href="/feedback/", className="ml-2", external_link=True)),
                    ],
                    align="center",# no_gutters=True,
                ),
                href="/",
            ),
            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            #dbc.Collapse(search_bar, id="navbar-collapse", navbar=True, is_open=False),
        ],
    color="dark", dark=True,
    ),

    dbc.Row(
        [
        dbc.Col(
            html.Div(
                #html.H1("Welcome to PlantApp"),
                #html.P('Dash converts Python classes into HTML'),
                dcc.Markdown('''
                    ## Welcome to PlantApp

                    The genomic resources for plants are rapidly expending. There are
                    hundreds of sequenced genomes across the genetic diversity of the plant kindgdom and within the 
                    genetic diversity of individual species. Each sequenced genome comes with useful human- and machine-
                    friendly gene annotation files. PlantApp uses [SQNce](https://github.com/eporetsky/SQNce), 
                    a SQLite-based database framework for parsing commonly used gene-annotation files. SQNce enables the generation 
                    of single, uniform, database that can is fast to generate and query. PlantApp aggregates publicly 
                    available data from sequenced plants and offers an easy-to-access and easy-to-download database. 

                    Both PlantApp and SQNce are currently under development and are constantly being updated. At this developmental
                    phase I am primarily relying on the [Phytozome](https://phytozome-next.jgi.doe.gov/) annotation files. 
                    For requests or issue report e-mail me at eporetsky at ucsd.edu or message me on Twitter @externelly.
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
    fluid=True  
)