from typing import Optional
from fastapi import FastAPI, Path, Query, HTTPException, Depends, status
from fastapi.responses import (
    StreamingResponse,
    Response,
    RedirectResponse,
    PlainTextResponse,
)

from fastapi.security import HTTPBasic, HTTPBasicCredentials

## Add Doc Security
from fastapi.openapi.docs import get_redoc_html, get_swagger_ui_html
from fastapi.openapi.utils import get_openapi

from enum import Enum
from pydantic import BaseModel
import IsoSpecPy as iso
import sys
import secrets
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, Draw
from rdkit.Chem.Draw import rdMolDraw2D

import re
from io import BytesIO
from PIL import Image

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt, ExactMolWt, HeavyAtomMolWt
import IsoSpecPy as iso
from matplotlib import pyplot as plt


class ImageFormat(str, Enum):
    svg = "svg"
    png = "png"


description = """
API to convert from different chemical formats to other formats using RDKit

## SMILES

You can **convert SMILES**.

* to InChi and InChi key.
* to images in SVG or PNG format.
* to isotopic distribution using [IsoSpecPy](http://matteolacki.github.io/IsoSpec/)
"""
expectedPassword = "<<API-Password>>"
app = FastAPI(
    title="Chemical Translator API",
    description=description,
    version="0.0.1",
    contact={
        "name": "Torsten Schindler",
        # "url": "https://www.roche.com/contact",
        "email": "Torsten.Schindler@roche.com",
    },
    # license_info={
    #     "name": "GPL v3",
    #     "url": "https://www.gnu.org/licenses/gpl-3.0.en.html",
    # },
)


## Initial Basic Auth
security = HTTPBasic()
## For document
def get_current_username(credentials: HTTPBasicCredentials = Depends(security)):
    correct_username = secrets.compare_digest(credentials.username, "user")
    correct_password = secrets.compare_digest(credentials.password, "<<API-Password>>")
    if not (correct_username and correct_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
            headers={"WWW-Authenticate": "Basic"},
        )
    return credentials.username

## For actual API => to aviod Authentication header from GCP's bearer token
def weak_authentication(credentials : dict):
    weak_key = credentials.get('weak_authen','')
    if weak_key != "temporary_work":
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    return 'authen_ok'

## Create Auth Doc
@app.get("/openapi.json", include_in_schema=False)
async def openapi(username: str = Depends(get_current_username)):
    return get_openapi(title=app.title, version=app.version, routes=app.routes)

@app.route('/')
@app.get("/docs", include_in_schema=False)
async def get_swagger_documentation(username: str = Depends(get_current_username)):
    return get_swagger_ui_html(openapi_url="/openapi.json", title="docs")

@app.get("/redoc", include_in_schema=False)
async def get_redoc_documentation(username: str = Depends(get_current_username)):
    return get_redoc_html(openapi_url="/openapi.json", title="docs")



def smile_to_2D_mol(smiles):
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    return mol


def get_png(smiles, width=200, height=200):
    mol = smile_to_2D_mol(smiles)
    img_io = BytesIO()
    pil_img = Draw.MolToImage(mol, size=(width, height), kekulize=True)
    pil_img.save(img_io, "PNG", quality=100)
    img_io.seek(0)
    return img_io


def get_svg(smiles, width=200, height=200):
    mol = smile_to_2D_mol(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg.replace("svg:", "")
    # Cut off the first line with xml tag
    svg = "\n".join(svg.split("\n")[1:])
    # Replace :svg from xmlns otherwise the SVG is not shown in bokeh
    svg = svg.replace("xmlns:svg=", "xmlns=")
    # Replace end of line by space
    svg = svg.replace("\n", " ")
    # Substitute multiple spaces and tabs by a single space
    svg = re.sub("\s+", " ", svg).strip()
    return svg


def create_isotop_plot(masses, probs, figsize=(2, 2), dpi=72, digits=2):
    """
    Create isotope distribution plot

    See annotate function in https://github.com/matrixx567/MassSpectraPlot/blob/master/MassSpectraPlot/MassSpectraPlot.py
    """
    # plt.figure(figsize=figsize, dpi=dpi)
    # ax = plt.subplot(111)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, frameon=True)
    ax.vlines(x=masses, ymin=[0], ymax=probs, linewidth=2.0)
    ax.set_xlabel("Masses")
    # Hide y axis
    ax.get_yaxis().set_visible(False)
    # Hide spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    # Annotate peaks
    for mass, prob in zip(masses, probs):
        # plt.text(mass, prob, f"{mass:.4f}", rotation=90)
        plt.annotate(
            f"{mass:.{digits}f}",
            xy=(mass, prob),
            xycoords="data",
            xytext=(mass, prob + 0.01),
            textcoords="data",
            # rotation=90,
            size=12,
            horizontalalignment="center",
            verticalalignment="bottom",
            # arrowprops=dict(
            #     arrowstyle="-",
            #     color="#808080",
            #     linewidth=0.4,
            #     shrinkA=0.05,
            #     shrinkB=1,
            # ),
        )
    fig.tight_layout()
    img_io = BytesIO()
    plt.savefig(img_io, format="png")
    plt.close()
    img_io.seek(0)
    return img_io


# @app.get("/", include_in_schema=False)
# async def docs_redirect():
#     """Redirect access to start page to docs

#     Exclude this endpoint from the docs
#     """
#     return RedirectResponse(url="docs")


@app.get("/smiles-to-inchi/{password}/{smiles:path}")
async def get_inchi_from_smiles(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to InChi"),
):
    """Get InChi string and key from SMILES

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    try:
        smiles = smiles.strip()
        mol = Chem.MolFromSmiles(smiles)
        result = {
            "smiles": smiles,
            "InChi": Chem.MolToInchi(mol),
            "InChiKey": Chem.MolToInchiKey(mol),
        }
        return result
    except:
        raise HTTPException(
            status_code=404, detail=f"Conversion of SMILES to InChi failed: {smiles}"
        )


@app.get("/smiles-to-inchi-plain/{password}/{smiles:path}", response_class=PlainTextResponse)
async def get_plain_inchi_from_smiles(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to InChi"),
):
    """Get plain text InChi from SMILES

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )    
    try:
        result = Chem.MolToInchi(Chem.MolFromSmiles(smiles.strip()))
        return result
    except:
        return "Failed"


@app.get("/smiles-to-inchikey-plain/{password}/{smiles:path}", response_class=PlainTextResponse)
async def get_plain_inchikey_from_smiles(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to InChi"),
):
    """Get plain text InChi key from SMILES

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )    
    try:
        result = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles.strip()))
        return result
    except:
        return "Failed"


@app.get("/inchi-to-inchikey/{password}/{inchi:path}")
async def get_inchi_key_from_string(
    inchi,
    password: str = Path(..., title="Password needed")):
    """Convert InChi to InChi key and SMILES

    Example:

    The InCHI String: `InChI=1S/C11H11N/c1-2-3-6-10-7-4-5-8-11(10)9-12/h2-5,7-8H,6H2,1H3/b3-2-`
    returns the InCHI Key: IVOJPQZCVBYKKF-IHWYPQMZSA-N
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )    
    try:
        inchi = inchi.strip()
        mol = Chem.MolFromInchi(inchi)
        result = {
            "InChi": inchi,
            "InChiKey": Chem.MolToInchiKey(mol),
            "smiles": Chem.MolToSmiles(mol),
        }
        return result
    except:
        raise HTTPException(
            status_code=404, detail=f"Conversion of InChi failed: {inchi}"
        )


@app.get("/inchi-to-inchikey-plain/{password}/{inchi:path}", response_class=PlainTextResponse)
async def get_plain_text_inchi_key_from_inchi(
    password: str = Path(..., title="Password needed"),
    inchi: str = Path(
        ..., title="InChi", description="InChi string to convert to InChi key"
    ),
):
    """Get plain text InChi key from InChi string

    Example:

    `InChI=1S/C11H11N/c1-2-3-6-10-7-4-5-8-11(10)9-12/h2-5,7-8H,6H2,1H3/b3-2-`
    """
    
    try:
        result = Chem.MolToInchiKey(Chem.MolFromInchi(inchi.strip()))
        return result
    except:
        return "Failed"


@app.get("/smiles-to-isodist/{password}/{smiles:path}")
async def get_iso_from_smiles(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to isotopic distribution"),
    probability: Optional[float] = Query(
        None,
        title="Coverage P",
        description="Coverage probability of isotopic distribution",
    ),
):
    """Get a given coverage P of the isotopic distribution from SMILES. See [IsoSpec](http://matteolacki.github.io/IsoSpec/)

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    try:
        smiles = smiles.strip()
        mol = Chem.MolFromSmiles(smiles)
        formula = CalcMolFormula(mol, separateIsotopes=False, abbreviateHIsotopes=True)
        # set default coverage to P = 99.99% = 0.9999 = .9999
        prob_to_cover = probability if probability else 0.9999
        sp = iso.IsoTotalProb(formula=formula, prob_to_cover=prob_to_cover)
        result = {
            "smiles": smiles,
            "formula": formula,
            "average_mol_weight": MolWt(mol),
            "exact_mol_weight": ExactMolWt(mol),
            "heavy_atom_mol_weight": HeavyAtomMolWt(mol),
            "prob_to_cover": prob_to_cover,
            # "isotop_dist": [{"mass": mass, "prob": prob} for mass, prob in sp],
            "masses": list(sp.masses),
            "probs": list(sp.probs),
        }
        return result
    except:
        raise HTTPException(
            status_code=404,
            detail=f"Conversion of SMILES to isotopic distribution failed: {smiles}",
        )


@app.get("/smiles-to-isoplot/{password}/{smiles:path}")
async def get_isoplot_from_smiles(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to isotopic distribution"),
    probability: Optional[float] = Query(
        0.95,
        title="Coverage P",
        description="Coverage probability of isotopic distribution",
        ge=0.0,
        le=1.0,
    ),
    width: Optional[float] = Query(
        3.2,
        title="Figure width",
        description="Width of the isotopic distribution figure in inch",
    ),
    height: Optional[float] = Query(
        2.4,
        title="Figure height",
        description="Height of the isotopic distribution figure in inch",
    ),
    dpi: Optional[float] = Query(
        100,
        title="Resolution of the figure",
        description="Resolution of the figure in dots per inch",
    ),
    digits: Optional[int] = Query(
        2, title="Digits", description="Number of digits displayed at the mass peaks"
    ),
):
    """Get a given coverage probability of the isotopic distribution from SMILES

    Example SMILES: `COc1cc(c2cn(Cc3cn4cc(CNC[C@@]56C[C@H](C5)C6)ccc4n3)nn2)c7cncn7c1`

    If a 640 x 480 pixel image of the isotopic distribution is required use
    * 6.4 inch width x 100 dots per inch = 640 pixel
    * 4.8 inch height x 100 dots per inch = 480 pixel
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    try:
        smiles = smiles.strip()
        mol = Chem.MolFromSmiles(smiles)
        formula = CalcMolFormula(mol, separateIsotopes=False, abbreviateHIsotopes=True)
        # set default coverage to P = 99.99% = 0.9999 = .9999
        prob_to_cover = probability if probability else 0.9999
        sp = iso.IsoTotalProb(formula=formula, prob_to_cover=prob_to_cover)
        masses = list(sp.masses)
        probs = list(sp.probs)
        buf = create_isotop_plot(
            masses, probs, figsize=(width, height), dpi=dpi, digits=digits
        )
        return StreamingResponse(buf, media_type="image/png")
    except:
        raise HTTPException(
            status_code=404,
            detail=f"Conversion of SMILES to isotopic distribution plot failed: {smiles}",
        )


@app.get("/smiles-to-image/{password}/{smiles:path}")
async def smiles_to_image(
    password: str = Path(..., title="Password needed"),
    smiles: str = Path(..., title="SMILES to convert to image"),
    width: Optional[int] = Query(200, title="Width", description="Image width", ge=1),
    height: Optional[int] = Query(
        200, title="Height", description="Image height", ge=1
    ),
    imageformat: Optional[ImageFormat] = Query(
        ImageFormat.png, title="Format", description="Image format", alias="format"
    ),
):
    """Get Image from SMILES

    Example SMILES:

    `O=C1CC2S[C@H]3C[C@]2(C=C1Br)C1=C(N3)C(=O)c2[nH]cc3c2C1=NCC3`

    `CC(C)C12C(=O)N3C4C(C(C3(C(=O)N1C)SS2)O)(C5=CC=CC=C5N4)C67C(C89C(=O)N(C(C(=O)N8C6NC1=CC=CC=C71)(SSSS9)CO)C)O`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    
    try:
        if imageformat == imageformat.svg:
            text = get_svg(smiles, width, height)
            return Response(text, media_type="image/svg+xml;utf8")
        elif imageformat == imageformat.png:
            image_bytes = get_png(smiles, width, height)
            return StreamingResponse(image_bytes, media_type="image/png")
        else:
            raise HTTPException(
                status_code=404, detail=f"Unsupported image format: {imageformat}"
            )
    except:
        raise HTTPException(
            status_code=404, detail=f"Conversion of SMILES failed: {smiles}"
        )


@app.get("/svg/{password}/{width}/{height}/{smiles:path}")
async def smiles_to_svg(
    password: str = Path(..., title="Password needed"),
    width: int = Path(..., title="Width of the image", ge=1),
    height: int = Path(..., title="Height of the image", ge=1),
    smiles: str = Path(..., title="SMILES to convert to PNG"),
):
    """Get SVG from SMILES

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )    
    try:
        text = get_svg(smiles, width, height)
        return Response(text, media_type="image/svg+xml;utf8")
    except:
        raise HTTPException(
            status_code=404, detail=f"Conversion of SMILES failed: {smiles}"
        )


@app.get("/png/{password}/{width}/{height}/{smiles:path}")
async def smiles_to_png(
    password: str = Path(..., title="Password needed"),
    width: int = Path(..., title="Width of the image", ge=1),
    height: int = Path(..., title="Height of the image", ge=1),
    smiles: str = Path(..., title="SMILES to convert to PNG"),
):

    if password != expectedPassword:
        raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect email or password",
                headers={"WWW-Authenticate": "Basic"},
            )
    """Get PNG from SMILES

    Example SMILES: `N#Cc1ccccc1C/C=C\C`
    """
    try:
        image_bytes = get_png(smiles, width, height)
        return StreamingResponse(image_bytes, media_type="image/png")
    except:
        raise HTTPException(
            status_code=404, detail=f"Conversion of SMILES failed: {smiles}"
        )
