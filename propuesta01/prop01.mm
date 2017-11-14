<map version="freeplane 1.3.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node ID="ID_1723255651" CREATED="1283093380553" MODIFIED="1506002316526"><richcontent TYPE="NODE">

<html>
  <head>
    
  </head>
  <body>
    <p style="text-align: center">
      Propuesta 1
    </p>
    <p style="text-align: center">
      Locomocion humanoide, contacto,
    </p>
    <p style="text-align: center">
      impacto y dinamica centroidal
    </p>
  </body>
</html>
</richcontent>
<hook NAME="MapStyle">

<map_styles>
<stylenode LOCALIZED_TEXT="styles.root_node">
<stylenode LOCALIZED_TEXT="styles.predefined" POSITION="right">
<stylenode LOCALIZED_TEXT="default" MAX_WIDTH="600" COLOR="#000000" STYLE="as_parent">
<font NAME="SansSerif" SIZE="10" BOLD="false" ITALIC="false"/>
</stylenode>
<stylenode LOCALIZED_TEXT="defaultstyle.details"/>
<stylenode LOCALIZED_TEXT="defaultstyle.note"/>
<stylenode LOCALIZED_TEXT="defaultstyle.floating">
<edge STYLE="hide_edge"/>
<cloud COLOR="#f0f0f0" SHAPE="ROUND_RECT"/>
</stylenode>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.user-defined" POSITION="right">
<stylenode LOCALIZED_TEXT="styles.topic" COLOR="#18898b" STYLE="fork">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.subtopic" COLOR="#cc3300" STYLE="fork">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.subsubtopic" COLOR="#669900">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.important">
<icon BUILTIN="yes"/>
</stylenode>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.AutomaticLayout" POSITION="right">
<stylenode LOCALIZED_TEXT="AutomaticLayout.level.root" COLOR="#000000">
<font SIZE="18"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,1" COLOR="#0033ff">
<font SIZE="16"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,2" COLOR="#00b439">
<font SIZE="14"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,3" COLOR="#990000">
<font SIZE="12"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,4" COLOR="#111111">
<font SIZE="10"/>
</stylenode>
</stylenode>
</stylenode>
</map_styles>
</hook>
<hook NAME="AutomaticEdgeColor" COUNTER="5"/>
<node TEXT="Introduccion" POSITION="right" ID="ID_1124539920" CREATED="1509553752442" MODIFIED="1509553769714">
<edge COLOR="#00ffff"/>
<node TEXT="Dificultades de la simulacion de la locomoci&#xf3;n humanoide y su aplicaci&#xf3;n a la biomecanica" ID="ID_1398657070" CREATED="1509553769718" MODIFIED="1509553844565" MAX_WIDTH="300" MIN_WIDTH="300">
<node TEXT="Modelo de colision y contacto" ID="ID_1174144689" CREATED="1509554825300" MODIFIED="1509554851193">
<node TEXT="Basados en penalizacion" ID="ID_792872692" CREATED="1509554590436" MODIFIED="1509554664698" MAX_WIDTH="200" MIN_WIDTH="200">
<node TEXT="Modelos" ID="ID_579533063" CREATED="1509554632273" MODIFIED="1509554656112">
<node TEXT="Hertz" ID="ID_1609837208" CREATED="1509554642285" MODIFIED="1509554650787"/>
<node TEXT="Hunt and Crossley" ID="ID_820877449" CREATED="1509554675696" MODIFIED="1509554685808"/>
<node TEXT="Cundall" ID="ID_616968783" CREATED="1509554668618" MODIFIED="1509554672060"/>
<node TEXT="Azad-Feathrestone" ID="ID_1642175480" CREATED="1509554714221" MODIFIED="1509554722818"/>
</node>
<node TEXT="ventajas y desventajas" ID="ID_91439937" CREATED="1509553863294" MODIFIED="1509564688583" MAX_WIDTH="200" MIN_WIDTH="200">
<node TEXT="dt pque&#xf1;os para una buena precision" ID="ID_1487982041" CREATED="1509553930264" MODIFIED="1509554089594" MAX_WIDTH="200" MIN_WIDTH="200"/>
<node TEXT="paremetrizacion del material" ID="ID_1196940644" CREATED="1509553949675" MODIFIED="1509554048138" MAX_WIDTH="200" MIN_WIDTH="200"/>
<node TEXT="penetraci&#xf3;n de cuerpos" ID="ID_1720165210" CREATED="1509553967824" MODIFIED="1509554048140" MAX_WIDTH="200" MIN_WIDTH="200"/>
<node TEXT="mecanismos para la conservacion de la energia" ID="ID_966274324" CREATED="1509553975372" MODIFIED="1509554048141" MAX_WIDTH="200" MIN_WIDTH="200"/>
<node TEXT="metodos adaptativos" ID="ID_1380918750" CREATED="1509554094478" MODIFIED="1509554297829"/>
<node TEXT="facil implementacion" ID="ID_372546053" CREATED="1509564918592" MODIFIED="1509641164857"/>
</node>
</node>
<node TEXT="Basados en conservacion del momento" ID="ID_17234522" CREATED="1509564717974" MODIFIED="1509564848668" MAX_WIDTH="200" MIN_WIDTH="200"/>
</node>
<node TEXT="Solucion de los modelos MBD en 3D" ID="ID_1079666304" CREATED="1509554866539" MODIFIED="1509554930718">
<node TEXT="ODE" ID="ID_1845662549" CREATED="1509554901086" MODIFIED="1509554903201">
<node TEXT="Sin articulaciones con singularidad que requieren de quaterniones, pueden usar cualquier metodo adaptativo y para sistemas rigidos" ID="ID_1084544026" CREATED="1509554941379" MODIFIED="1509555193950" MAX_WIDTH="200" MIN_WIDTH="200"/>
</node>
<node TEXT="DAE" ID="ID_1692842040" CREATED="1509554903656" MODIFIED="1509554911573">
<node TEXT="Requieren el uso de quaterniones y es necesario el uso de esquemas especiales, frecuentemente Euler implicito" ID="ID_524937713" CREATED="1509555200670" MODIFIED="1509555261652" MAX_WIDTH="200" MIN_WIDTH="200"/>
</node>
</node>
</node>
</node>
<node TEXT="Simulacion" POSITION="right" ID="ID_1543090780" CREATED="1506002323768" MODIFIED="1506002335241">
<edge COLOR="#ff0000"/>
<node TEXT="Motores fisico" ID="ID_192242353" CREATED="1506002335250" MODIFIED="1506002341447">
<node TEXT="ODE" ID="ID_1650791484" CREATED="1506002341452" MODIFIED="1506002346309"/>
<node TEXT="Featherstone&apos;s approach" ID="ID_894574550" CREATED="1506002347130" MODIFIED="1506002363387"/>
</node>
<node TEXT="6D-vectors" ID="ID_319523912" CREATED="1506002366521" MODIFIED="1506002373266"/>
</node>
<node TEXT="Modelo de contacto" POSITION="right" ID="ID_583940179" CREATED="1506002415886" MODIFIED="1506002423918">
<edge COLOR="#0000ff"/>
<node TEXT="Hertz" ID="ID_539195732" CREATED="1506002423924" MODIFIED="1506002427862"/>
<node TEXT="Mozah" ID="ID_1257937700" CREATED="1506002428536" MODIFIED="1506002444649"/>
</node>
<node TEXT="Dinamica centroidal" POSITION="right" ID="ID_1772502658" CREATED="1506002447392" MODIFIED="1506002460102">
<edge COLOR="#00ff00"/>
<node TEXT="Momento lineal" ID="ID_565957192" CREATED="1506002471041" MODIFIED="1506002479256"/>
<node TEXT="Momento angular" ID="ID_1411318346" CREATED="1506002460105" MODIFIED="1506002469971"/>
</node>
<node TEXT="Generaci&#xf3;n de trayectorias" POSITION="right" ID="ID_1011506672" CREATED="1506002489529" MODIFIED="1506002495915">
<edge COLOR="#ff00ff"/>
</node>
</node>
</map>
